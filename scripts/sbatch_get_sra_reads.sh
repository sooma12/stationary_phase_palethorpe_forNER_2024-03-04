#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=get_sra_reads
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --array=1-6%7
#SBATCH --ntasks=6
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=soo.m@northeastern.edu

# sbatch_get_sra_reads.sh

module load anaconda3/2021.11
source activate BINF-12-2021

ID_LIST=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/data/rawreads/srr_id.list

id=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_LIST |  awk '{print $1}')

fasterq-dump --split-3 $id -O /work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/data/rawreads/
