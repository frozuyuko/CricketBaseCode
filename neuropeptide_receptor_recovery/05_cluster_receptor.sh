set -euo pipefail

7TM=~/CricketBaseCode/genomeannotation/data/interproscan/gpcr_7tm.pep.fa
NO7TM="~/CricketBaseCode/genomeannotation/data/interproscan/gpcr_no7tm.pep.fa"            
OUT="~/CricketBaseCode/genomeannotation/data/interproscan/gpcr_cluster"
mkdir -p "$OUT"

pool="$OUT/pool.fa"
#combine 7TM and non-7TM candidates into one pool
cat "$7TM" ${NO7TM:+"$NO7TM"} > "$pool"

# cluster: -c identity, -aS coverage of the shorter seq, -d 0 keep full names
cd-hit -i "$pool" -o "$OUT/reps.fa" -c 0.95 -aS 0.5 -d 0 -T 0 -M 0

echo "representatives -> $OUT/reps.fa  (centroids)"
echo "cluster membership -> $OUT/reps.fa.clstr"