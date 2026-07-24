set -euo pipefail

BLAST_DB="~/CricketBaseCode/genomeannotation/data/blastdb"
blastp -query ~/CricketBaseCode/genomeannotation/data/gbi_knownneuropeptide.fa \
  -db "$BLAST_DB" \
  -evalue 1e-5 \
  -num_threads 8 \
  -out ~/CricketBaseCode/genomeannotation/data/neuropeptide_receptor_blastp.txt

blastp -query ~/CricketBaseCode/genomeannotation/data/newneuropeptide_neurostresspep.fa \
-db "$BLAST_DB" \
  -evalue 1e-5 \
  -num_threads 8 \
  -out ~/CricketBaseCode/genomeannotation/data/neuropeptide_receptor_blastp.txt


# Then manually check output and for those that are not found, use rna-seq guided recovery 






