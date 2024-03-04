# Reanalyzing Palethorpe et al (JBac 2022) BfmRS stationary phase data for sRNAs
## MWS started March 4th 2024

Palethorpe S, Farrow JM, Wells G, Milton ME, Actis LA, Cavanagh J, Pesci EC. 2022.Acinetobacter baumannii Regulates Its Stress Responses via the BfmRS Two-Component Regulatory System. J Bacteriol 204:e00494-21.https://doi.org/10.1128/jb.00494-21

### Get data
Relevant lines from methods section:
"The processed reads were mapped to the strain ATCC 17961 draft genome sequence (http://www.patricbrc.org; genome ID 470.2202 [70])"
"Raw Illumina sequencing reads were deposited in the NCBI Sequence Read Archive under BioProject ID PRJNA780533, with accession numbers SRX13141184, SRX13141185, and SRX13141186 for strain ATCC 17961 and SRX13141187, SRX13141188, and SRX13141189 for strain 17961-ΔbfmR."
"Gene names and locus tags in the text have been updated to the identifiers found in the completed genome sequence for strain ATCC 17961 (GenBank accession numbers CP065432 [chromosome], CP065433 [pAB17961-1], and CP065434 [pAB17961-2]) (72)."

```bash
# Get compute node for downloading data
bash ~/request_compute_node.bash 

cd /work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04

mkdir -p genome
init_agent # custom function for dealing with git ssh
touch .gitignore
echo '/genome' >.gitignore
git cmp '' # custom function: add-commit-push

# Download reference genome
cat genome_extensions.list 
#fna
#features.tab
#gff
#fnn
#frn

for ext in `cat genome_extensions.list`; do wget -a genome/wget.log -N "ftp://ftp.bvbrc.org/genomes/470.2202/470.2202.$ext"; done
# fna downloaded correctly; all others yielded 404.  Tried individually and still got 404.

mv 470* genome

# Download RNA-seq data
mkdir -p data/rawreads
module load anaconda3/2021.11
source activate BINF-12-2021
touch data/rawreads/srr_id.list
cat data/rawreads/srr_id.list
#SRX13141184
#SRX13141185
#SRX13141186
#SRX13141187
#SRX13141188
#SRX13141189
# "Raw Illumina sequencing reads were deposited in the NCBI Sequence Read Archive under BioProject ID PRJNA780533, with accession numbers SRX13141184, SRX13141185, and SRX13141186 for strain ATCC 17961 and SRX13141187, SRX13141188, and SRX13141189 for strain 17961-ΔbfmR."
# Sbatch this command...
sbatch scripts/sbatch_get_sra_reads.sh 

```

### Fastqc to find adapter sequence

```bash
cd ./scripts
sbatch sbatch_2_fastqc.sh
```
