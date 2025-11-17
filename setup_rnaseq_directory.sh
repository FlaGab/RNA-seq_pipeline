#!/bin/bash
# This script sets up the directory for RNA-seq analysis by creating necessary folders.
# main folder is defined in settings.sh which is placed in the same directory as this script.
source ./settings.txt

# Create directory structure"
mkdir -p "$main_folder"/{raw_data,QC_raw_data,aligned_data,counts,differential_expression,plots,reports}
if [ "$trimming" = true ]; then
    mkdir -p "$main_folder/QC_trimmed_data"
    mkdir -p "$main_folder/QC_trimmed_data"
fi 

if [ "$pe" = true ]; then
    mkdir -p "$main_folder"/trimmed_data/{paired,unpaired}
fi

echo "Directory structure created under $main_folder"

# place fq.gz files in raw_data folder
echo "Please place your .fq.gz files in the 'raw_data' folders located at: $main_folder/raw_data"
echo "if sequencing is paired-end, ensure that files are named with _1 and _2 suffixes for forward and reverse reads respectively."