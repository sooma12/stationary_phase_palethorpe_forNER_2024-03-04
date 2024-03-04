#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=fastqc
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.output
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.error
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

# based on 2024 Jan RNA-seq script from my ∆pbpG/∆lpsB

echo "Starting fastqc SBATCH script $(date)"

echo "Loading environment and tools"
#fastqc requires OpenJDK/19.0.1
module load OpenJDK/19.0.1
module load fastqc/0.11.9

FASTQDIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/data/rawreads
OUT_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/data/fastqc_output
SCRIPT_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/scripts
LOG_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs

mkdir -p $FASTQDIR $OUT_DIR $SCRIPT_DIR

echo "Running fastqc in directory $FASTQDIR"
fastqc $FASTQDIR/*.fastq

echo "Cleaning up logs and output files"
mv $SCRIPT_DIR/fastqc_rnaseq* $LOG_DIR
mkdir -p $OUT_DIR/fastqc_html $OUT_DIR/fastqc_zip
mv $FASTQDIR/*fastqc.html $OUT_DIR/fastqc_html
mv $FASTQDIR/*fastqc.zip $OUT_DIR/fastqc_zip
