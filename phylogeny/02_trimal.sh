#!/usr/bin/env bash
# trimAl v1.5.rev1

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p alignment

echo "[1/3] main"
trimal -in alignment/comb_aln.fasta -out alignment/comb_trim.fasta -gt 0.2 -cons 40

echo "[2/3] npr"
trimal -in alignment/npr_aln.fasta -out alignment/npr_trim.fasta -gt 0.2 -cons 40

echo "[3/3] ilp"
trimal -in alignment/ilp_aln.fasta -out alignment/ilp_trim.fasta -gt 0.2 -cons 40