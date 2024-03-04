#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=trimmomatic
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

# sbatch sbatch_3_trim.sh

# Load trimmomatic and java
module load trimmomatic/0.39
module load oracle_java/jdk1.8.0_181

# Directories for inputs and outputs.  NOTE: MUST END WITH / FOR SHELL SUBSTITUTION BELOW TO WORK
BASE_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04
FASTQ_INDIR=${BASE_DIR}/data/rawreads/
PAIRED_OUTDIR=${BASE_DIR}/data/trimmed/paired/
UNPAIRED_OUTDIR=${BASE_DIR}/data/trimmed/unpaired/

# path to NU Discovery cluster's Trimmomatic program folder with Illumina adapters
PATH_TO_TRIMMOMATIC="/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2"

## Initialize variables to contain file suffixes and output paths
leftInSuffix="_1.fastq"
rightInSuffix="_2.fastq"
leftOutSuffix="_trimmed_R1.fastq"
rightOutSuffix="_trimmed_R2.fastq"

mkdir -p $PAIRED_OUTDIR $UNPAIRED_OUTDIR

# Loop through all the left-read fastq files in INDIR
for leftInFile in ${FASTQ_INDIR}*${leftInSuffix}
do
  # Remove path from filename
  pathRemoved="${leftInFile/"$FASTQ_INDIR"/}"
  # Remove left-read suffix, leaving the name of the sample
  sampleName="${pathRemoved/$leftInSuffix/}"
  echo Trimming $sampleName
  # Test with echos; comment this out before final use
#  echo "$FASTQ_INDIR$sampleName$leftInSuffix"
#  echo "$FASTQ_INDIR$sampleName$rightInSuffix"
#  echo "$PAIRED_OUTDIR$sampleName$leftOutSuffix"
#  echo "$UNPAIRED_OUTDIR$sampleName$leftOutSuffix"
#  echo "$PAIRED_OUTDIR$sampleName$rightOutSuffix"
#  echo "$UNPAIRED_OUTDIR$sampleName$rightOutSuffix"
  # Use sample name derived from shell replacement to trim left AND right reads
  java -jar /shared/centos7/trimmomatic/0.39/trimmomatic-0.39.jar PE \
  -threads 1 -phred33 \
  $FASTQ_INDIR$sampleName$leftInSuffix \
  $FASTQ_INDIR$sampleName$rightInSuffix \
  $PAIRED_OUTDIR$sampleName$leftOutSuffix \
  $UNPAIRED_OUTDIR$sampleName$leftOutSuffix \
  $PAIRED_OUTDIR$sampleName$rightOutSuffix \
  $UNPAIRED_OUTDIR$sampleName$rightOutSuffix \
  HEADCROP:0 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36 \
  ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10
done