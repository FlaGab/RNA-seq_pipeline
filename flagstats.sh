#!/bin/bash
# This script runs samtools stats on all .bam files generated during alignment step
source ./settings.txt

#find all .bam files in the aligned data directory
files=$(ls "${main_folder}/aligned_data"/*.bam)

#loop through each file
for file in $files
do
  #extract the sample name from the file name
  basename=$(echo $file | cut -d '.' -f 1)
  #run samtools flagstat on the file and save the output to a new file
  samtools stats $file > "${basename}_stats.txt"
done