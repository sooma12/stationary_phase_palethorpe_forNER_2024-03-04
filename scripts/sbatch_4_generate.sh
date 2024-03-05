#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=genomeGenerate
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch sbatch_4_generate.sh

module load star/2.7.11a

# Variables
GENOME_REF_DIR_OUT=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref
FASTA_IN=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/470_2202/470.2202.fna
GTF_IN=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/470_2202/470_2202.gtf

# STAR requires the output directory be pre-made
mkdir -p $GENOME_REF_DIR_OUT

# STAR time
# Program recommended `--genomeSAindexNbases 9` after running with default value 14

STAR --runMode genomeGenerate \
--genomeDir /work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/ \
--genomeFastaFiles $FASTA_IN \
--sjdbGTFfile $GTF_IN \
--genomeSAindexNbases 9 \
--sjdbGTFfeatureExon transcript
