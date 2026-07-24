# CricketBaseCode

A comprehensive sequence- and structure-based resource for the neuropeptide-receptor interactome of a non-model insect, *Gryllus bimaculatus*

Febrina Margaretha, Mika Sakamoto, Hitomi Seike, Shinji Nagata, Yasukazu Nakamura, and Takako Mochizuki

Code accompanying the genome re-annotation, neuropeptide and GPCR curation, structure-based interaction screen, and receptor phylogeny underlying [CricketBase](https://cricket.annotation.jp).

## Pipeline overview

```
genomeannotation/            re-annotate the genome, extract GPCR candidates
        |
neuropeptide_receptor_recovery/   recover missing/partial neuropeptides and receptors
        |
structureprediction/         fold every peptide-receptor pair, score each interface
        |
phylogeny/                   place candidate receptors against known ligand-specificity references
```

### 1. `genomeannotation/`

Structural and functional re-annotation of the chromosome-scale assembly.

- `01_braker4_annotate.sh`: BRAKER4 (ETP mode: GeneMark-ETP + AUGUSTUS) gene prediction, using sixteen paired-end RNA-seq libraries (7 developmental-stage SRR + 9 dissected-tissue DRR) and Arthropoda protein evidence from OrthoDB v12. --> preferably run in slurm
- `02_interproscan6.sh`: InterProScan 6 (Nextflow workflow) functional domain annotation of the predicted proteome.
- `03_busco610.sh`: BUSCO v6.1.0 completeness against the insecta_odb12 and arthropoda_odb12 lineages.
- `samples.csv`, `config.ini`: BRAKER4 sample sheet and run configuration.

### 2. `neuropeptide_receptor_recovery/`

Recovery of neuropeptide precursors and GPCR candidates missed or left incomplete by the initial annotation.

- `01_neuropeptideblast.sh`: BLASTP of known and NeuroStresPep-derived neuropeptide precursors against the predicted proteome.
- `02_neuropeptidestringtie.sh`: tblastn plus StringTie/TransDecoder-guided recovery of precursors not found by BLASTP.
- `03_deeptmhmm.sh`: DeepTMHMM validation of seven-transmembrane (7TM) topology for rhodopsin- and secretin-class receptor candidates.
- `04_stringtie_sequence_recovery.sh`: StringTie/TransDecoder rescue of receptor loci with an incomplete 7TM topology.
- `05_cluster_receptor.sh`: CD-HIT clustering (95% identity, 50% coverage of the shorter sequence) of 7TM candidates to non-redundant representatives.

### 3. `structureprediction/`

Structure prediction of every peptide-receptor pair with AlphaFold3 and Boltz-2, and interface confidence scoring.

- `01_fastatojson.py`: build AlphaFold3 input JSON from precursor/receptor FASTA.
- `02_submitaf3msa_{neuropeptide,rhodopsin,secretin}.sh`: AlphaFold3 MSA generation (SLURM array jobs).
- `03_af3msatoboltz2yaml.py`: convert AlphaFold3 MSAs to Boltz-2 paired-input YAML (`af3tools`).
- `04_af3run_{rhodopsin,secretin}.sh`: AlphaFold3 v3.0.1 structure prediction (SLURM array jobs, GPU).
- `05_boltz2run.sh`: Boltz-2 structure prediction (SLURM array job, GPU).
- `06_structuralscoring.sh`: interface scoring driver: `ipsae_ptm.py` (ipSAE), `af3_plddt_ptm.py` / `boltz2_plddt_ptm.py` (interface pLDDT, size-weighted and balanced), all at a PAE cutoff of 10 A and a distance cutoff of 10 A.
- `createptm.sh`: prepares post-translational-modification-aware structure inputs (pyroglutamate, C-terminal amidation, disulfide bonds, N-glycosylation) ahead of Boltz-2 prediction.
- `fix_pairedreceptorinput.py`: validates and corrects PTM positions, disulfide bonds, and dimer glycosylation/linkage in the paired AlphaFold3 input JSONs.
- `peptide_ptms.tsv`: per-peptide PTM annotation table consumed by the scripts above.

### 4. `phylogeny/`

Maximum-likelihood placement of candidate receptors against reference receptors of known ligand specificity. Builds three trees: **main** (240 seqs, 75 reference + 165 cricket 7TM representatives), **npr** (76 reference-only neuropeptide GPCRs), and **ilp** (6 seqs, cricket insulin-like receptor with insect/human InsR/IGF1R).

- `01_mafft.sh`: MAFFT v7.526 (`--auto --reorder`) alignment of all three sequence sets (`input/combined_165_ref.fasta`, `input/neuropeptide_receptors.fasta`, `input/ilp_tree_input.fasta`).
- `02_trimal.sh`: trimAl v1.5.rev1 (`-gt 0.2 -cons 40`) trimming of all three alignments.
- `03_iqtree3.sh`: IQ-TREE 3.1.2 maximum-likelihood trees, 1000 ultrafast bootstrap and 1000 SH-aLRT replicates on each: main (`LG+F+G4`), npr (`MFP` ModelFinder, best-fit `LG+F+I+R6`), ilp (`LG+G4`).
- `confidentmap.py`: prunes the main tree to its 200-tip representatives-only view for the browser, and derives the phylogeny-confident receptor set (SH-aLRT >= 80 and UFBoot >= 95).

## Resource

The annotation, curated neuropeptide and receptor sets, ranked predicted complexes, and receptor phylogeny produced by this pipeline are served through [CricketBase](https://cricket.annotation.jp).

## Citation

If you use this code or the associated CricketBase resource, please cite:

Margaretha F, Sakamoto M, Seike H, Nagata S, Nakamura Y, Mochizuki T. A comprehensive sequence- and structure-based resource for the neuropeptide-receptor interactome of the non-model insect *Gryllus bimaculatus*. (in preparation).

```bibtex
@article{margaretha_cricketbase,
  author  = {Margaretha, Febrina and Sakamoto, Mika and Seike, Hitomi and Nagata, Shinji and Nakamura, Yasukazu and Mochizuki, Takako},
  title   = {A comprehensive sequence- and structure-based resource for the neuropeptide-receptor interactome of the non-model insect {G}ryllus bimaculatus},
  note    = {in preparation},
  url     = {https://cricket.annotation.jp}
}
```
