#!/bin/bash
#SBATCH --job-name=af3_msa_rho
#SBATCH -p l40s
#SBATCH -A l40s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:0
#SBATCH --output=/home/febrinam-pg/finalcricketnrp/logs/msa/af3_msa_sec_%A_%a.out
#SBATCH --error=/home/febrinam-pg/finalcricketnrp/logs/msa/af3_msa_sec_%A_%a.err
#SBATCH --array=1-68

set -euo pipefail

IMAGE_PATH="/lustre12/software/alphafold3/alphafold3-v3.0.1.sif"
ALPHAFOLD3DIR="/app/alphafold"
MODEL_DIR="/home/febrinam-pg/tools/alphafold3/models"

INPUT_DIR="/home/febrinam-pg/finalcricketnrp/json/receptor/secretin"
MSA_OUT_BASE="/home/febrinam-pg/finalcricketnrp/af3/msa/receptor/secretin"

PEPTIDE_LIST="${INPUT_DIR}/secretinreceptor_list.txt"

# Get the JSON file name for this array index
JSON_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PEPTIDE_LIST}")
if [ -z "${JSON_FILE}" ]; then
  echo "No JSON_FILE for task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 1
fi

# Full input JSON path
INPUT_JSON_PATH="${INPUT_DIR}/${JSON_FILE}"

# Peptide name without .json suffix, e.g. Gbim.scaffold_6G0000170.1_input
PEPTIDE_BASENAME="${JSON_FILE%.json}"

# Output dir for this peptide (AF3 will create a subdir named after \"name\" field)
OUTPUT_DIR="${MSA_OUT_BASE}/${PEPTIDE_BASENAME}"
mkdir -p "${OUTPUT_DIR}"

echo "Task ${SLURM_ARRAY_TASK_ID}: ${INPUT_JSON_PATH} -> ${OUTPUT_DIR}"

apptainer exec "${IMAGE_PATH}" \
  python3.11 "${ALPHAFOLD3DIR}/run_alphafold.py" \
    --model_dir="${MODEL_DIR}" \
    --json_path="${INPUT_JSON_PATH}" \
    --output_dir="${OUTPUT_DIR}" \
    --norun_inference
