#!/bin/bash
#SBATCH --job-name=af3pred_sec
#SBATCH -p l40s
#SBATCH -A l40s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --output=/home/febrinam-pg/finalcricketnrp/logs/predict/af3pred_sec_%A_%a.out
#SBATCH --error=/home/febrinam-pg/finalcricketnrp/logs/predict/af3pred_sec_%A_%a.err
#SBATCH --array=1-119

set -euo pipefail
IMAGE_PATH="/lustre12/software/alphafold3/alphafold3-v3.0.1.sif"
ALPHAFOLD3DIR="/app/alphafold"
MODEL_DIR="/home/febrinam-pg/tools/alphafold3/models"

CLASS="secretin"
EXPECTED=68
PAIRS_BASE="/home/febrinam-pg/finalcricketnrp/af3/pairedreceptorinput/${CLASS}"
OUT_BASE="/home/febrinam-pg/finalcricketnrp/af3/structureoutput/${CLASS}"

mapfile -t PEPTIDES < <(ls -d "${PAIRS_BASE}"/*/ | xargs -n1 basename | sort)
PEPTIDE="${PEPTIDES[$((SLURM_ARRAY_TASK_ID-1))]:-}"
if [ -z "${PEPTIDE}" ]; then echo "no peptide folder for task ${SLURM_ARRAY_TASK_ID}"; exit 0; fi

INPUT_DIR="${PAIRS_BASE}/${PEPTIDE}"
OUTPUT_DIR="${OUT_BASE}/${PEPTIDE}"
mkdir -p "${OUTPUT_DIR}"
N=$(ls -1 "${INPUT_DIR}"/*.json | wc -l)
if [ "${N}" -lt "${EXPECTED}" ]; then echo "WARNING: ${PEPTIDE} has ${N}, expected ${EXPECTED}"; fi
echo "Task ${SLURM_ARRAY_TASK_ID}: ${PEPTIDE}, folding ${N} pairings -> ${OUTPUT_DIR}"

apptainer exec --nv "${IMAGE_PATH}" \
  python3.11 "${ALPHAFOLD3DIR}/run_alphafold.py" \
    --model_dir="${MODEL_DIR}" \
    --input_dir="${INPUT_DIR}" \
    --output_dir="${OUTPUT_DIR}" \
    --norun_data_pipeline