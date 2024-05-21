#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=featureCounts
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

echo "Loading environment and tools"
module load subread/2.0.6

source ./config.cfg
echo "featureCounts file name: $COUNTS_FILE found in $COUNTS_OUTDIR"
echo "genome GTF reference file: $GENOME_GTF"
echo ".bam input files were found in: $STAR_MAPPED_DIR"

mkdir -p $COUNTS_OUTDIR

# Run featureCounts on all BAM files from STAR
# -t flag specifies features in column 3 of GTF to count; default is exon.
featureCounts \
-a $GENOME_GTF \
-o $COUNTS_OUTDIR/$COUNTS_FILE \
-p \
-t transcript \
$STAR_MAPPED_DIR/*.bam
