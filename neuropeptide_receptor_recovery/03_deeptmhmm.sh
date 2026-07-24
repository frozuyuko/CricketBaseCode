set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate newdeeptmhmm
python -c "import torch; print(f'torch {torch.__version__}  cuda_avail={torch.cuda.is_available()}')"

TOOL_DIR=~/tools/DeepTMHMM-Academic-License-v1.0
IPR_DIR=~/CricketBaseCode/genomeannotation/data/interproscan
RHOD_OUT=$IPR_DIR/deeptmhmm_rhodopsin
SECR_OUT=$IPR_DIR/deeptmhmm_secretin

cd "$TOOL_DIR"
# clean-rerun pattern (your working-script style): wipe then let predict.py recreate
rm -rf "$RHOD_OUT" "$SECR_OUT"

python predict.py --fasta "$IPR_DIR/rhodopsin.fa" --output-dir "$RHOD_OUT"
python predict.py --fasta "$IPR_DIR/secretin.fa"  --output-dir "$SECR_OUT"

for d in "$RHOD_OUT" "$SECR_OUT"; do
  echo "== TM helix tally: $d =="
  awk 'BEGIN{RS=">"} NR>1{split($0,a,"\n"); t=a[3]; n=gsub(/M+/,"x",t); print a[1]"\t"n" TM"}' \
    "$d/predicted_topologies.3line"
done
