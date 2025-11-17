#!/bin/bash
#this script create the reference genome index for STAR aligner
#load variables from settings file
source ./settings.txt

# compute index for STAR
echo "Indexing genome with STAR..."
STAR --runThreadN $threads \
     --runMode genomeGenerate \
     --genomeDir "indexed_reference_genome" \
     --genomeFastaFiles "$path_to_genome_sequence" \
     --sjdbGTFfile "$path_to_genome_annotation" \
     --sjdbOverhang ${sjdbOverhang} \
     --genomeSAindexNbases ${genomeSAindexNbases}
     