#!/usr/bin/env bash
# IQ-TREE 3.1.2
#
set -euo pipefail
cd "$(dirname "$0")/.."

THREADS=8
mkdir -p trees

echo "[1/3] main - LG+F+G4"
iqtree -s alignment/comb_trim.fasta -m LG+F+G4 -B 1000 -alrt 1000 \
    -T AUTO -ntmax $THREADS --prefix trees/comb -redo

echo "[2/3] npr - MFP (ModelFinder; best-fit was LG+F+I+R6)"
iqtree -s alignment/npr_trim.fasta -m MFP -B 1000 -alrt 1000 \
    -T AUTO -ntmax $THREADS --prefix trees/npr -redo

echo "[3/3] ilp - LG+G4, 4 threads"
iqtree -s alignment/ilp_trim.fasta -m LG+G4 -B 1000 -alrt 1000 \
    -T AUTO -ntmax 4 --prefix trees/ilptree -redo
cp trees/ilptree.treefile trees/ilp_receptor_tree.treefile

echo "Done. To regenerate the tree and the"
echo "confident receptor set, run python confidentmap.py."