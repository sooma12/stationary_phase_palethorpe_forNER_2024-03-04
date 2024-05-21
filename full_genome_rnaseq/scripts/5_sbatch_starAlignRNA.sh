#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=alignRNA
#SBATCH --time=04:00:00
#SBATCH --array=1-9%10
#SBATCH --ntasks=9
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/logs/%x-%j-%a.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/logs/%x-%j-%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch 5_sbatch_starAlignRNA.sh
echo "Loading tools"
module load star/2.7.11a
source ./config.cfg

echo "Reference genome directory: $GENOME_REF_DIR"
echo "sample sheet located at $SAMPLE_SHEET_PATH"
echo "alignment output directory containing .bam files: ${DATA_DIR}/mapped/"

# STAR requires the output directory be pre-made
mkdir -p ${DATA_DIR}/mapped

# Array job!  Used my sample sheet technique from 2023-11-13 breseq for this.
# NOTE!!! sample sheet prep is moved to 4_prep_sample_sheet_for_starAlign.sh.
# Run that script before this one!

# sed and awk read through the sample sheet and grab each whitespace-separated value
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET_PATH |  awk '{print $1}')
r1=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET_PATH |  awk '{print $2}')
r2=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET_PATH |  awk '{print $3}')

echo "Running STAR on files $r1 and $r2"

# STAR time
# note alignIntronMax=1 to disallow introns for bacteria
# BAM sorting seems to be memory intensive.  Set an available amount of memory using limitBAMsortRAM 
STAR --runMode alignReads \
--genomeDir $GENOME_REF_DIR \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $r1 $r2 \
--runThreadN 4 \
--alignIntronMax 1 \
--limitBAMsortRAM 5000000000 \
--outFileNamePrefix ${DATA_DIR}/mapped/${name}
