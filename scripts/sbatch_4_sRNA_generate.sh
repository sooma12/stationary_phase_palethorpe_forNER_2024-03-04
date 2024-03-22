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

## Usage: sbatch sbatch_4_sRNA_generate.sh

module load star/2.7.11a

# Variables
GENOME_REF_DIR_OUT=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/CP065432_1_srnas
FASTA_IN=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/CP065432_1_srnas/CP065432.1.fasta
GTF_IN=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/ref/CP065432_1_srnas/srnas_17961_24-03-22.gtf

# STAR requires the output directory be pre-made
mkdir -p $GENOME_REF_DIR_OUT

# STAR time
# Program recommended `--genomeSAindexNbases 9` after running with default value 14

STAR --runMode genomeGenerate \
--genomeDir $GENOME_REF_DIR_OUT \
--genomeFastaFiles $FASTA_IN \
--sjdbGTFfile $GTF_IN \
--genomeSAindexNbases 9 \
--sjdbGTFfeatureExon sRNA
