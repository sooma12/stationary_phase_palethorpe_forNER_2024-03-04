#!/bin/bash
# prep_sample_sheet_for_starAlign.sh
# Makes a sample_sheet.txt containing sample ID and R1 and R2 filepaths
# Assumes each sample file is in the format: WT_1_S84_trimmed_R2.fastq
# Script uses a split on "_" to grab the sample ID, e.g. WT_1.  Must modify this if sample file names are different!
# Example output line:  WT_1 /path/to/WT_1_R1.fastq /path/to/WT_1_R2.fastq
# Usage: bash script_5_prep_sample_sheet_for_STAR.sh <path/to/fastqs> <path/to/samplesheet>

SAMPLE_SHEET=${2}

# Create .list files with R1 and R2 fastqs.  Sort will put them in same orders, assuming files are paired
find ${1} -maxdepth 1 -name "*.fastq" | grep -e "_R1" | sort > R1.list
find ${1} -maxdepth 1 -name "*.fastq" | grep -e "_R2" | sort > R2.list

if [ -f "${SAMPLE_SHEET}" ] ; then
  rm "${SAMPLE_SHEET}"
fi

# make sample sheet from R1 and R2 files.  Format on each line looks like (space separated):
# WT_1 /path/to/WT_1_R1.fastq /path/to/WT_1_R2.fastq
# from sample sheet, we can access individual items from each line with e.g. `awk '{print $3}' sample_sheet.txt`

paste R1.list R2.list | while read R1 R2 ;
do
    outdir_root=$(basename ${R2} | cut -f1 -d"_")
    sample_line="${outdir_root} ${R1} ${R2}"
    echo "${sample_line}" >> $SAMPLE_SHEET
done

rm R1.list R2.list
