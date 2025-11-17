#!/bin/bash
# This scripts execute alignment of RNA-seq reads to a reference genome using STAR.
source ./settings.txt

#check if the reference genome index exists, if not create it
if [ ! -d "$path_to_indexed_reference_genome" ] || [ -z "$(ls -A $path_to_indexed_reference_genome)" ] || [ $path_to_indexed_reference_genome = "false" ]; then
    echo "Reference genome index not found. Creating index..."
    bash ./genome_index.sh
    path_to_indexed_reference_genome="indexed_reference_genome"
fi

if [ "$pe" != true ]; then
    # List of filenames
    # Get all unique base names
    # Sort and use 'uniq' to ensure each base name is processed once
    BASENAMES=$(ls "${main_folder}"/trimmed_data/*.fq.gz | \
        xargs -n 1 basename | \
        sed 's/fq.gz//' | sort | uniq)

    # Loop over each filename and run STAR alignment
    for filename in "${BASENAMES[@]}"; do
        echo "Aligning $filename..."
        STAR \
        --runThreadN $threads \
        --readFilesCommand gunzip -c \
        --genomeDir  $path_to_indexed_reference_genome\
        --readFilesIn "${main_folder}/trimmed_data/${filename}.fq.gz" \
        --outFileNamePrefix "${main_folder}/aligned_data/${filename}" \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outTmpDir "${main_folder}/tmp" 
    done
else
    # List of filenames
    # Get all unique base names (e.g., name1, name2)
    # Sort and use 'uniq' to ensure each base name is processed once
    BASENAMES=$(ls -1 "${main_folder}"/trimmed_data/paired/*_1P.fq.gz 2>/dev/null | \
        xargs -n 1 basename | \
        sed 's/_1P.fq.gz//' | sort | uniq)

    # Loop over each filename and run STAR alignment
    for filename in "${BASENAMES[@]}"; do
        echo "Aligning $filename..."
        STAR \
        --runThreadN $threads \
        --readFilesCommand gunzip -c \
        --genomeDir  $path_to_indexed_reference_genome\
        --readFilesIn "${main_folder}/trimmed_data/paired/${filename}_1P.fq.gz" "${main_folder}/trimmed_data/paired/${filename}_2P.fq.gz" \
        --outFileNamePrefix "${main_folder}/aligned_data/${filename}" \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outTmpDir "${main_folder}/tmp"
    done
fi  

echo "All files processed."