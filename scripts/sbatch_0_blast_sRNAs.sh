#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=blast_sRNAs
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.output
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.error
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

# Files/dirs
WORK_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04
SRNA_FA=${WORK_DIR}/17978-srnas.fasta
DB_17961=${WORK_DIR}/ref/CP065432_ref


# Load conda environment with BLAST
# shellcheck disable=SC2046
module load ncbi-blast+/2.13.0

blastn -query $SRNA_FA -db ${DB_17961}/CP065432.1.fasta \
    -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 4
