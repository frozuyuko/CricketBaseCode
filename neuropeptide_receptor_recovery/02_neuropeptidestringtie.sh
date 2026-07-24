set -euo pipefail 

MISSING_PEP = ~/CricketBaseCode/genomeannotation/data/missingneuropeptide.fa
GENOME = ~/CricketBaseCode/genomeannotation/data/genome.fa.masked
BAMS = ~/CricketBaseCode/genomeannotation/data/bam/*.sorted.bam
PAD = 2000
MIN_ORF = 50
OUT = ~/CricketBaseCode/genomeannotation/data/recover_missing_out
THREADS = 8
mkdir -p "$OUT"

tblastn -query "$MISSING_PEP" \
  -db "$BLAST_DB" \
  -evalue 1e-5 \
  -num_threads 8 \
  -out "$OUT/missingneuropeptide.txt"

# Determine loci 

# Miniprot + Hisat2 + StringTie + TransDecoder to recover missing neuropeptides

for t in makeblastdb tblastn stringtie miniprot gffread bedtools TransDecoder.LongOrfs TransDecoder.Predict python3; do
  command -v "$t" >/dev/null 2>&1 || { echo "[error] missing tool: $t"; exit 1; }
done
[[ -s "$GENOME" ]]      || { echo "[error] genome not found: $GENOME"; exit 1; }
[[ -s "$MISSING_PEP" ]] || { echo "[error] missing-precursor FASTA not found: $MISSING_PEP"; exit 1; }

makeblastdb -in "$GENOME" -dbtype nucl -out "$OUT/genome_db" >/dev/null

tblastn -query "$MISSING_PEP" -db "$OUT/genome_db" -evalue 1e-5 \
  -num_threads "$THREADS" -outfmt 6 > "$OUT/tblastn.tsv"
echo "tblastn hits: $(wc -l < "$OUT/tblastn.tsv")"

python3 - "$OUT/tblastn.tsv" "$PAD" "$OUT/best_hit.bed" <<'PY'
import sys, collections
tsv, pad, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
best = {}
for ln in open(tsv):
    c = ln.rstrip("\n").split("\t")
    q, s, sstart, send, bit = c[0], c[1], int(c[8]), int(c[9]), float(c[11])
    lo, hi = min(sstart, send), max(sstart, send)
    if q not in best or bit > best[q][3]:
        best[q] = (s, lo, hi, bit)
with open(out, "w") as o:
    for q, (s, lo, hi, bit) in sorted(best.items()):
        start = max(0, lo - 1 - pad)
        end = hi + pad
        o.write(f"{s}\t{start}\t{end}\t{q}\n")
print(f"loci (best hit per query, padded {pad} bp): {len(best)}")
PY

sort -k1,1 -k2,2n "$OUT/best_hit.bed" > "$OUT/best_hit.sorted.bed"
bedtools merge -i "$OUT/best_hit.sorted.bed" -c 4 -o distinct > "$OUT/loci.bed"
echo "merged loci: $(wc -l < "$OUT/loci.bed")"

gtfs=()
IFS=',' read -ra bamlist <<< "$BAMS"
for bam in "${bamlist[@]}"; do
  g="$OUT/$(basename "$bam" .bam).gtf"
  [[ -s "$g" ]] || stringtie "$bam" -p "$THREADS" -o "$g"
  gtfs+=("$g")
done
printf '%s\n' "${gtfs[@]}" > "$OUT/gtf_list.txt"
[[ -s "$OUT/stringtie_merged.gtf" ]] || stringtie --merge -o "$OUT/stringtie_merged.gtf" "$OUT/gtf_list.txt"

python3 - "$OUT/loci.bed" "$OUT/stringtie_merged.gtf" "$OUT/overlap.gtf" <<'PY'
import sys, re
loci_bed, merged_gtf, out_gtf = sys.argv[1], sys.argv[2], sys.argv[3]
iv = {}
for ln in open(loci_bed):
    if not ln.strip():
        continue
    c = ln.split("\t")
    iv.setdefault(c[0], []).append((int(c[1]), int(c[2])))
for k in iv:
    iv[k].sort()

def overlaps(chrom, s, e):
    return any(not (e < a or s > b) for a, b in iv.get(chrom, []))

tid_re = re.compile(r'transcript_id "([^"]+)"')
keep = set()
for ln in open(merged_gtf):
    if ln.startswith("#"):
        continue
    p = ln.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] != "transcript":
        continue
    if overlaps(p[0], int(p[3]), int(p[4])):
        m = tid_re.search(p[8])
        if m:
            keep.add(m.group(1))

with open(merged_gtf) as f, open(out_gtf, "w") as o:
    for ln in f:
        if ln.startswith("#"):
            o.write(ln)
            continue
        m = tid_re.search(ln)
        if m and m.group(1) in keep:
            o.write(ln)
print(f"transcripts overlapping missing-precursor loci: {len(keep)}")
PY

miniprot -t "$THREADS" --gff "$GENOME" "$MISSING_PEP" > "$OUT/miniprot.gff3"

gffread -w "$OUT/transcripts.fa" -g "$GENOME" "$OUT/overlap.gtf"
NTX=$(grep -c '^>' "$OUT/transcripts.fa" || echo 0)
echo "candidate transcripts: $NTX"
[[ "$NTX" -gt 0 ]] || { echo "[error] no candidate transcripts at the missing-precursor loci"; exit 1; }

TransDecoder.LongOrfs -t "$OUT/transcripts.fa" -m "$MIN_ORF" --output_dir "$OUT/td"

makeblastdb -in "$MISSING_PEP" -dbtype prot -out "$OUT/missing_pep_db" >/dev/null
blastp -query "$OUT/td/longest_orfs.pep" -db "$OUT/missing_pep_db" \
  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads "$THREADS" \
  > "$OUT/td/longest_orfs.blastp.tsv"

TransDecoder.Predict -t "$OUT/transcripts.fa" --single_best_only --no_refine_starts \
  --retain_blastp_hits "$OUT/td/longest_orfs.blastp.tsv" --output_dir "$OUT/td"
