set -euo pipefail

PROJ="~/CricketBaseCode"; TOOLS="${HOME}/tools"
RESULTS="~/CricketBaseCode/genomeannotation/output"
BRAKER_GTF="$RESULTS/braker.gtf"
LOCI_BED="~/CricketBaseCode/genomeannotation/data/interproscan/gpcr_no7tm.loci.bed"
GENOME="~/CricketBaseCode/genomeannotation/data/genome.fa.masked"
BAM_DIR="~/CricketBaseCode/genomeannotation/data/bam"

WORK="~/CricketBaseCode/genomeannotation/data/stringtie"
ST_PS="${WORK}/stringtie_per_sample"; ST_MG="${WORK}/stringtie_merged"
RESCUE_DIR="${WORK}/rescue"; TD_DIR="${WORK}/transdecoder"; PEP_DIR="${WORK}/peptides"
MERGED_GTF="${ST_MG}/merged.gtf"; ERROR_BED="${RESCUE_DIR}/errorgenes.bed"
OVERLAP_GTF="${RESCUE_DIR}/merged.overlap.gtf"
RESCUE_FA="${TD_DIR}/rescue_transcripts.fa"; TD_WORK="${TD_DIR}/td_workdir"
RAW_PEP="${TD_WORK}/rescue_transcripts.fa.transdecoder.pep"
FINAL_PEP="${PEP_DIR}/rescue.pep.fa"

ACCESSIONS=(SRR14026720 SRR14026721 SRR14026722 SRR14026723 SRR14026724 SRR14026725 SRR14026726
            DRR358356 DRR358357 DRR358358 DRR358359 DRR358360 DRR358361 DRR358362 DRR358363 DRR358364)

mkdir -p "${ST_PS}" "${ST_MG}" "${RESCUE_DIR}" "${TD_DIR}" "${PEP_DIR}" "${WORK}/logs"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

source "~/miniconda3/etc/profile.d/conda.sh"
[[ -d "${TOOLS}/TransDecoder-TransDecoder-v6.0.0" ]] && \
  export PATH="${TOOLS}/TransDecoder-TransDecoder-v6.0.0:${TOOLS}/TransDecoder-TransDecoder-v6.0.0/util:${PATH}"
[[ -d "${TOOLS}/gffread-0.12.7.Linux_x86_64" ]] && \
  export PATH="${TOOLS}/gffread-0.12.7.Linux_x86_64:${PATH}"

echo "Node: $(hostname)  Start: $(date)  Threads: ${THREADS}"
stringtie --version; gffread --version

for f in "${BRAKER_GTF}" "${LOCI_BED}" "${GENOME}"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing: $f"; exit 1; }
done
echo "Rescue loci: $(wc -l < ${LOCI_BED})"
for ACC in "${ACCESSIONS[@]}"; do
  [[ -s "${BAM_DIR}/${ACC}.sorted.bam" ]] || { echo "[ERROR] missing BAM: ${ACC}"; exit 1; }
done
echo "All 16 reused BAMs present."
cp -f "${LOCI_BED}" "${ERROR_BED}"

echo ">>> [Phase 1] per-sample StringTie (reused BAMs, BRAKER guide)"
PER=()
for ACC in "${ACCESSIONS[@]}"; do
  ST="${ST_PS}/${ACC}.stringtie.gtf"
  if [[ -s "${ST}" ]]; then echo "  [skip] ${ACC}"
  else echo "  stringtie ${ACC}"
       stringtie "${BAM_DIR}/${ACC}.sorted.bam" -G "${BRAKER_GTF}" -o "${ST}" -p "${THREADS}"; fi
  PER+=("${ST}")
done

echo ">>> [Phase 2] StringTie --merge"
LIST="${ST_MG}/gtf_list.txt"; printf "%s\n" "${PER[@]}" > "${LIST}"
[[ -s "${MERGED_GTF}" ]] || stringtie --merge -G "${BRAKER_GTF}" -o "${MERGED_GTF}" -p "${THREADS}" "${LIST}"
echo "  Merged transcripts: $(awk '$3=="transcript"' ${MERGED_GTF} | wc -l)"

echo ">>> [Phase 3] subset merged GTF to the 48 loci"
if [[ ! -s "${OVERLAP_GTF}" ]]; then
ERROR_BED="${ERROR_BED}" MERGED_GTF="${MERGED_GTF}" OVERLAP_GTF="${OVERLAP_GTF}" python3 - <<'PYEOF'
import os, re
iv = {}
for ln in open(os.environ["ERROR_BED"]):
    if not ln.strip(): continue
    c = ln.split("\t"); iv.setdefault(c[0], []).append((int(c[1]), int(c[2])))
for k in iv: iv[k].sort()
def ov(ch, s, e):
    return any(not (e < a or s > b) for a, b in iv.get(ch, []))
tid = re.compile(r'transcript_id "([^"]+)"'); keep = set()
for ln in open(os.environ["MERGED_GTF"]):
    if ln.startswith("#"): continue
    p = ln.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] != "transcript": continue
    if ov(p[0], int(p[3]), int(p[4])):
        m = tid.search(p[8])
        if m: keep.add(m.group(1))
with open(os.environ["MERGED_GTF"]) as f, open(os.environ["OVERLAP_GTF"], "w") as o:
    for ln in f:
        if ln.startswith("#"): o.write(ln); continue
        m = tid.search(ln)
        if m and m.group(1) in keep: o.write(ln)
print(f"Overlapping transcripts: {len(keep)}")
PYEOF
fi
echo "  Overlap transcripts: $(awk '$3=="transcript"' ${OVERLAP_GTF} | wc -l)"

echo ">>> [Phase 4] gffread + TransDecoder"
[[ -s "${RESCUE_FA}" ]] || gffread -w "${RESCUE_FA}" -g "${GENOME}" "${OVERLAP_GTF}"
NTX=$(grep -c '^>' "${RESCUE_FA}" || echo 0); echo "  Rescue transcripts: ${NTX}"
[[ "${NTX}" -gt 0 ]] || { echo "[ERROR] no rescue transcripts produced"; exit 1; }

cd "${TD_DIR}"
if [[ ! -s "${RAW_PEP}" ]]; then
  TransDecoder.LongOrfs -t "${RESCUE_FA}" -O "${TD_WORK}"
  TransDecoder.Predict   -t "${RESCUE_FA}" --single_best_only --no_refine_starts -O "${TD_WORK}"
fi
for ext in pep bed cds gff3; do
  ln -sf "td_workdir/rescue_transcripts.fa.transdecoder.${ext}" \
         "${TD_DIR}/rescue_transcripts.fa.transdecoder.${ext}"
done

if [[ ! -s "${FINAL_PEP}" ]]; then
RAW_PEP="${RAW_PEP}" FINAL_PEP="${FINAL_PEP}" python3 - <<'PYEOF'
import os, re
raw = os.environ["RAW_PEP"]; out = os.environ["FINAL_PEP"]; n = 0
with open(raw) as f, open(out, "w") as o:
    cur = None; buf = []
    def flush():
        global n
        if cur is None: return
        s = re.sub(r"[\*\s]", "", "".join(buf))
        if s: o.write(f">{cur}\n{s}\n"); n += 1
    for ln in f:
        if ln.startswith(">"): flush(); cur = ln[1:].split()[0]; buf = []
        else: buf.append(ln.strip())
    flush()
print(f"Sanitized sequences: {n}")
PYEOF
fi