#!/bin/bash
#SBATCH --job-name=boltz_str
#SBATCH -p b200
#SBATCH -A b200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --gres=gpu:1
#SBATCH --nodelist=b200-04
#SBATCH --output=/home/febrinam-pg/finalcricketnrp/logs/boltz2str/str_%A_%a.out
#SBATCH --error=/home/febrinam-pg/finalcricketnrp/logs/boltz2str/str_%A_%a.err

set -euo pipefail
ROOT=/home/febrinam-pg/finalcricketnrp
CLASS="${CLASS:-rhodopsin}"
mkdir -p "${ROOT}/logs/boltz2str"

source /home/febrinam-pg/miniconda3/etc/profile.d/conda.sh
conda activate boltz_l40s_fresh
NVRTC_DIR="/home/febrinam-pg/miniconda3/envs/boltz_l40s_fresh/lib/python3.11/site-packages/nvidia/cu13/lib"
export LD_LIBRARY_PATH="${NVRTC_DIR}:${LD_LIBRARY_PATH:-}"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2

YAML_ROOT="${ROOT}/boltz2/pairedreceptorinput/${CLASS}"
OUT_ROOT="${ROOT}/boltz2/structureoutput/${CLASS}"

mapfile -t PEPS < <(ls -d "${YAML_ROOT}"/*/ 2>/dev/null | xargs -n1 basename | sort)
PEP="${PEPS[$((SLURM_ARRAY_TASK_ID-1))]:-}"
if [ -z "${PEP}" ]; then
  echo "no peptide folder for task ${SLURM_ARRAY_TASK_ID} (class ${CLASS}); exiting cleanly"
  exit 0
fi

PEP_DIR="${YAML_ROOT}/${PEP}"
PEP_OUT_DIR="${OUT_ROOT}/${PEP}"
mkdir -p "${PEP_OUT_DIR}"
echo "Task ${SLURM_ARRAY_TASK_ID}  class=${CLASS}  peptide=${PEP}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "PEP_DIR=${PEP_DIR}"
echo "PEP_OUT_DIR=${PEP_OUT_DIR}"

cd "${PEP_DIR}"
shopt -s nullglob
n_ok=0; n_fail=0
for cfg in *.yaml; do
  [ -f "${cfg}" ] || continue
  echo "[${CLASS} ${PEP}] boltz predict ${cfg}"
  if boltz predict \
        --model boltz2 \
        --accelerator gpu \
        --devices 1 \
        --num_workers 20 \
        --override \
        --out_dir "${PEP_OUT_DIR}" \
        "${cfg}"; then
    n_ok=$((n_ok+1))
  else
    echo "FAILED ${CLASS}/${PEP}/${cfg}" >&2
    n_fail=$((n_fail+1))
  fi
done

echo "${CLASS}/${PEP} done: ${n_ok} ok, ${n_fail} failed."