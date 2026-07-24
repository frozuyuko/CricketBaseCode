#!/usr/bin/env bash
# MAFFT v7.526
#
# Aligns the sequence sets for all three CricketBase phylogeny trees:
#   main - 240 seqs (75 reference + 165 cricket 7TM representatives)
#   npr  - 76 reference-only neuropeptide GPCRs
#   ilp  - 6 seqs (cricket insulin-like receptor + insect/human InsR/IGF1R)

set -euo pipefail
cd "$(dirname "$0")/.."

THREADS=8
mkdir -p alignment trees

echo "[1/3] main - 240 seqs"
mafft --auto --thread $THREADS --reorder input/combined_165_ref.fasta \
    > alignment/comb_aln.fasta 2> trees/comb_mafft.log

echo "[2/3] npr - 76 seqs"
mafft --auto --thread $THREADS --reorder input/neuropeptide_receptors.fasta \
    > alignment/npr_aln.fasta 2> trees/npr_mafft.log

echo "[3/3] ilp - 6 seqs"
mafft --auto --thread 4 --reorder input/ilp_tree_input.fasta \
    > alignment/ilp_aln.fasta 2> trees/ilp_mafft.log