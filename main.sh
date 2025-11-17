#!/bin/bash
# This script execute autonomously all the steps from raw sequencing data to count matrix generation
# It assumes that the necessary tools (FastQC, Trimmomatic, STAR, featureCounts) are installed and available in the PATH
# It also assumes that the reference genome and annotation files are available
# Please modify and save the parameters in setting.txt before running the script

# 1 - initialize analysis directory
./setup_analysis_directory.sh

# 2 - quality check of raw data
./QC_raw.sh

# 3 - trimming of raw data
./trimming.sh

# 4 - quality check of trimmed data
./QC_trimmed.sh

# 5 - alignment to reference genome
./alignment.sh

# 6 - post-alignment statistics
./flagstats.sh

# 7 - generation of count matrix
./featureCounts.sh