#!/bin/bash
# 4_prep_sample_sheet_for_starAlign.sh
# Makes a sample_sheet.txt containing sample ID and R1 and R2 filepaths
# Assumes each sample file is in the format: SRR16949322_trimmed_R1.fastq
# Example output line:  SRR16949322 /path/to/SRR16949322_R1.fastq /path/to/SRR16949322_R2.fastq
# Usage: bash script_5_prep_sample_sheet_for_STAR.sh <path/to/fastqs> <path/to/samplesheet>

SAMPLE_SHEET=${2}

# Create .list files with R1 and R2 fastqs.  Sort will put them in same orders, assuming files are paired
find ${1} -maxdepth 1 -name "*.fastq" | grep -e "_R1" | sort > R1.list
find ${1} -maxdepth 1 -name "*.fastq" | grep -e "_R2" | sort > R2.list

if [ -f "${SAMPLE_SHEET}" ] ; then
  rm "${SAMPLE_SHEET}"
fi

# make sample sheet from R1 and R2 files.
# from sample sheet, we can access individual items from each line with e.g. `awk '{print $3}' sample_sheet.txt`

paste R1.list R2.list | while read R1 R2 ;
do
    outdir_root=$(basename ${R2} | cut -f1 -d"_")
    sample_line="${outdir_root} ${R1} ${R2}"
    echo "${sample_line}" >> $SAMPLE_SHEET
done

rm R1.list R2.list
