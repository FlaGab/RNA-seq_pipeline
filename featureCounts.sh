#!/bin/bash
#this script performs read counting using featureCounts
source ./settings.txt

#find all .bam files in aligned folder
files=$(ls "${main_folder}/aligned_data"/*.bam)
# detect sequencing type
if [ "$pe" = true ]; then
    featureCounts -p --countReadPairs -t exon -g gene_id -a "${path_to_genome_annotation}" -o "${main_folder}/counts/fC_count_matrix.txt" "${main_folder}/aligned_data"/*.bam
else
    featureCounts -t exon -g gene_id -a "${path_to_genome_annotation}" -o "${main_folder}/counts/fC_count_matrix.txt" "${main_folder}/aligned_data"/*.bam
fi