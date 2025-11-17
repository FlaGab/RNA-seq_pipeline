#!/bin/bash
# This script execute quality cjeck on raw RNA-seq data using FastQC.
source ./settings.txt

#check if 0_QC_raw_data folder exists, if not create it
if [ ! -d "$main_folder/QC_raw_data" ]; then
  mkdir "$main_folder/QC_raw_data"
fi 

# Run FastQC on all .fq.gz files in raw_data folder and output to QC_raw_data folder
fastqc -o "$main_folder/QC_raw_data" "$main_folder/raw_data/"*.fq.gz -t ${threads}

# Run multiQC on the FastQC output
if [ ! -d "$main_folder/QC_raw_data/multiQC" ]; then
  mkdir "$main_folder/QC_raw_data/multiQC"
fi
multiqc -o "$main_folder/QC_raw_data/multiQC" "$main_folder/QC_raw_data"

echo "FastQC and MultiQC scan completed. Results are in $main_folder/QC_raw_data"