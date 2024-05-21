#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=fastqc_rnaseq
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/full_genome_rnaseq/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/full_genome_rnaseq/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

echo "Starting fastqc SBATCH script $(date)"

echo "Loading environment and tools"
#fastqc requires OpenJDK/19.0.1
module load OpenJDK/19.0.1
module load fastqc/0.11.9
source ./config.cfg

echo "fastq file directory: $FASTQDIR"
echo "output directory for fastQC: $FASTQC_OUT_DIR"
echo "script directory: $SCRIPT_DIR"

mkdir -p $FASTQDIR $FASTQC_OUT_DIR $SCRIPT_DIR

echo "Running fastqc in directory $FASTQDIR"
fastqc $FASTQDIR/*.fastq

echo "Cleaning up logs and output files"
mkdir -p $SCRIPT_DIR/logs
mv $SCRIPT_DIR/fastqc_rnaseq* $SCRIPT_DIR/logs
mkdir -p $FASTQC_OUT_DIR/fastqc_html $FASTQC_OUT_DIR/fastqc_zip
mv $FASTQDIR/*fastqc.html $FASTQC_OUT_DIR/fastqc_html
mv $FASTQDIR/*fastqc.zip $FASTQC_OUT_DIR/fastqc_zip
