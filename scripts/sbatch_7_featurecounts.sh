#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=featureCounts
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

# Usage: sbatch 6_sbatch_featurecounts.sh
BASE_DIR=/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04
OUT_DIR=${BASE_DIR}/data/featurecounts
# GTF file for sRNAs in 17961
GTF_REF=${BASE_DIR}/ref/CP065432_1_srnas/srnas_17961_24-03-22.gtf
STAR_DIR=${BASE_DIR}/data/mapped/

echo "Loading environment and tools"
module load subread/2.0.6

mkdir -p $OUT_DIR

# Run featureCounts on all BAM files from STAR
# For detecting sRNAs, need -t flag (indicates what feature in GTF column 3; default is exon.)
featureCounts \
-a $GTF_REF \
-o $OUT_DIR/srna_counts.txt \
-p \
-t sRNA \
$STAR_DIR/*.bam
