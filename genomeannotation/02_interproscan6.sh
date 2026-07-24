source ~/miniconda3/etc/profile.d/conda.sh
conda activate braker4
export NXF_VER=25.04.6
export NXF_SINGULARITY_CACHEDIR=~/cricket/braker4pipeline/.singularity_cache

RESULTS=~/CricketBaseCode/genomeannotation/output/braker4
OUTDIR=~/CricketBaseCode/genomeannotation/output

mkdir -p $OUTDIR
cd $OUTDIR

# Please clean sequence with stop codons '*' 
INPUT_FA=~/CricketBaseCode/genomeannotation/data/Gbi_protein.fa

# Run InterProScan 6 
nextflow run ~/tools/interproscan6/main.nf \
    --input $INPUT_FA \
    --datadir $OUTDIR \
    --outdir $OUTDIR \
    --outprefix Gbi_protein \
    --formats "tsv,gff3,json,xml" \
    --goterms \
    --pathways \
    --cpus 96 \
    --max-workers 10 \
    -profile singularity \
    -resume