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

