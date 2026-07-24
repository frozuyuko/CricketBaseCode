#!/bin/bash
#SBATCH --job-name=boltz_prep
#SBATCH -p b200
#SBATCH -A b200
#SBATCH --nodelist=b200-04
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --output=/home/febrinam-pg/finalcricketnrp/logs/boltz/prep_%j.out
#SBATCH --error=/home/febrinam-pg/finalcricketnrp/logs/boltz/prep_%j.err

set -euo pipefail
ROOT=/home/febrinam-pg/finalcricketnrp
SC="${ROOT}/script"
mkdir -p "${ROOT}/logs/boltz"

echo "[0] reproduce a3m from AF3 JSONs (peptide_af3 + receptor _input dirs)"
python3 "${SC}/reproduce_boltz_msa.py" \
  --pep_af3_dir  "${ROOT}/json/peptide_af3" \
  --pep_src_dir  "${ROOT}/json/peptide" \
  --rec_af3_base "${ROOT}/af3/msa/receptor" \
  --out          "${ROOT}/boltz2/msa" \
  --apply

echo "[1] extract PTM table (2-Cys disulfides auto-snap; flags in this log)"
python3 "${SC}/extract_ptms.py" \
  --in_dir "${ROOT}/json/peptide" \
  --out    "${SC}/peptide_ptms.tsv" \
  --apply

echo "[2] monomer YAMLs (peptides chain B + PTMs, dimers B+C, receptors chain A)"
python3 "${SC}/make_boltz_monomer_yamls.py" \
  --pep_msa_dir  "${ROOT}/boltz2/msa/peptide" \
  --rec_msa_base "${ROOT}/boltz2/msa/receptor" \
  --pep_json_dir "${ROOT}/json/peptide" \
  --out_base     "${ROOT}/boltz2/yamlmonomer" \
  --ptms         "${SC}/peptide_ptms.tsv" \
  --apply

echo "PREP DONE"