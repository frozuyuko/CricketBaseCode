# Please install busco 6.1.0.
set -euo pipefail

PROTEOME="${PROTEOME:?BRAKER4 proteome FASTA (braker.aa)}"
OUT="${OUT:-busco_out}"
THREADS="${THREADS:-16}"

for lineage in insecta_odb12 arthropoda_odb12; do
  busco -i "$PROTEOME" -m proteins -l "$lineage" \
    -o "${lineage}" --out_path "$OUT" -c "$THREADS" -f
done

echo "BUSCO summaries -> $OUT/{insecta_odb12,arthropoda_odb12}/short_summary.*.txt"