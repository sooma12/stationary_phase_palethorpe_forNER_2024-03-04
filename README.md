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

mkdir -p ref
init_agent # custom function for dealing with git ssh
touch .gitignore
echo '/ref' >.gitignore
git cmp '' # custom function: add-commit-push

# Download reference genome

#cat genome_extensions.list 
##fna
##features.tab
##gff
##fnn
##frn
#for ext in `cat genome_extensions.list`; do wget -a genome/wget.log -N "ftp://ftp.bvbrc.org/genomes/470.2202/470.2202.$ext"; done
# fna downloaded correctly; all others yielded 404.  Tried individually and still got 404.
# IMPORTANT NOTE: I discovered that the RefSeq annotations are not available. This necessitated use of the PATRIC annotations...
##  Fixed version
cat genome_extensions.list 
#fna
#PATRIC.fna
#PATRIC.features.tab
#PATRIC.gff

for ext in `cat ref/genome_extensions.list`; do wget -a ref/wget.log -N "ftp://ftp.bvbrc.org/genomes/470.2202/470.2202.$ext"; done

mkdir -p ref/470_2202/
mv 470* ref/470_2202/

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

The authors reported that Illumina TruSeq was used for library preparation.

```bash
cd ./scripts
sbatch sbatch_2_fastqc.sh
```

Fastqc reports "Illumina universal adapter sequence".  I guessed which adapter file (that comes with trimmomatic) and trimmed two samples:

```bash
module load trimmomatic/0.39
module load oracle_java/jdk1.8.0_181
java -jar /shared/centos7/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 rawreads/SRR16949318_1.fastq rawreads/SRR16949318_2.fastq trimmed/paired/SRR16949318_1_paired.fastq trimmed/unpaired/SRR16949318_1_unpaired.fastq trimmed/paired/SRR16949318_2_paired.fastq trimmed/unpaired/SRR16949318_2_unpaired.fastq HEADCROP:0 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36 ILLUMINACLIP:/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10
```

Then, I ran fastqc: 
```bash
module load OpenJDK/19.0.1
module load fastqc/0.11.9
fastqc ./*.fastq
```

The resulting fastqc html files showed that adapters were removed, indicating that `TruSeq3-PE-2.fa` is the adapter file to use.

`sbatch_3_trim.sh` was the run to trim the all fastq files.

### Genome building - FIRST ATTEMPT; SCRATCH THIS

bv-brc.org seems not to have a gtf or gff file for the RefSeq fna file.  So, I'm going to use the fna along with the .PATRIC.gff.

Convert the gff file to gtf using gffread:
```bash
# Note, I set up gffread as a module in 2024-01_rnaseq_pbpGlpsB
module load gffread
REF_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/470_2202
gffread  $REF_DIR/470.2202.PATRIC.gff -T -o $REF_DIR/470_2202.gtf
```

Run STAR aligner in mode genomegenerate using sbatch_4_generate.sh.
Mar 5, 2024: This script threw the following error: 
`Fatal INPUT FILE error, no valid exon lines in the GTF file: /work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/470_2202/470_2202.gtf
Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.
`

Checked GTF file and the fna file to se if chromosome names agreed; they differed by 'accn|' in the GTF.  Edited gtf file with `sed`:
`sed 's/accn|//g' /work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/470_2202/470_2202.gtf >/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/470_2202/470_2202_chrrename.gtf`

### Genome building v2 - sRNAs only

Spoke to EG and decided on a different approach - begin with the complete genome (i.e. the full circular chromosome).
Find sRNAs by BLAST
Make an annotation file with just this.

1. Get sRNA list in FASTA format using `kroger_srnas_to_fasta.py`
2. Download 17961 fasta to `./ref/CP065432` from GenBank via `wget -O CP065432.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP065432.1&rettype=fasta"`
3. Blast sRNA fasta list against the CP065432.1 fasta file
- `conda activate /work/geisingerlab/conda_env/blast_corr`  I had installed BLAST previously
- Make BLAST database in `./ref/CP065432`: `makeblastdb -in CP065432.1.fasta -out CP065432 -parse_seqids -dbtype nucl`
- Run `sbatch scripts/sbatch_0_blast_sRNAs.sh`


### Alignment to reference

Prepare sample sheet using:
`bash scripts/script_5_prep_sample_sheet_for_STAR.sh /work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/data/trimmed/paired /work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/data/trimmed/paired/sample_sheet.txt`

Run STAR aligner using `sbatch_6_STARalign.sh`