#!/bin/bash
# Download files from
# Load conda environment with:
# `module load anaconda3`
# `conda activate /work/geisingerlab/conda_env/blast_corr`

module load anaconda3
source /shared/centos7/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate /work/geisingerlab/conda_env/blast_corr

source ./config.cfg

echo "Fetching fasta and gff3 files for accession numbers: " "${ACCESSION_ARRAY[@]}"
echo "Saving files to $REFERENCE_DIR"
echo "Genome GTF file output: $GENOME_GTF"

mkdir -p $REFERENCE_DIR

for i in "${ACCESSION_ARRAY[@]}"; do
  fasta_file=${i}.fasta
  gff3_file=${i}.gff3
  wget -O "$REFERENCE_DIR/$fasta_file" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${i}"
  wget -O "$REFERENCE_DIR/$gff3_file" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${i}"
done

for file in $REFERENCE_DIR/CP*.gff3; do
  cat "$file" > $REFERENCE_DIR/merged_17961.gff3
done

gffread  $REFERENCE_DIR/merged_17961.gff3 -T -o $GENOME_GTF
