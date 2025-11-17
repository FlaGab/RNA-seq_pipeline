#!/bin/bash
# This script executes quality check after trimming using FastQC.
source ./settings.txt

# Check if trimming is enabled
if [ "$trimming" = true ]; then
    echo "Starting quality check on trimmed data..."
    if [ "$pe" = true ]; then
        fastqc -t $threads -o "${main_folder}/QC_trimmed_data" "${main_folder}/trimmed_data/paired"/*.fq.gz
    else
        fastqc -t $threads -o "${main_folder}/QC_trimmed_data" "${main_folder}/trimmed_data/"/*.fq.gz
    fi
    multifastqc "${main_folder}/QC_trimmed_data" -o "${main_folder}/QC_trimmed_data/multiqc_report_trimmed"
    
    echo "Quality check on trimmed data completed. Reports are saved in ${main_folder}/QC_trimmed_data."
else
    echo "Trimming is disabled. Skipping quality check on trimmed data."
fi  