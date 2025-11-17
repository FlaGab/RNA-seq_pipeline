#!/bin/bash
# This script perform trimming of RNA-seq raw data using Trimmomatic.
source ./settings.txt

# Check if raw_data directory exists
if [ ! -d "${main_folder}/raw_data" ]; then
    echo "Error: raw_data directory does not exist in $main_folder. Please create it and add your .fq.gz files."
    exit 1
fi

# end this script if trimming is set to false
if [ "$trimming" = false ]; then
    echo "Trimming is set to false. Exiting trimming script."
    exit 0
fi

# set ILLUMINACLIP parameters
if [ -z "$adapters_file" ] || [ "$adapters_file" = false ]; then
    echo "No adapters file specified. Skipping ILLUMINACLIP step."
    ILLUMINACLIP=""
else
    ILLUMINACLIP="ILLUMINACLIP:${adapters_file}:${seed_mismatches}:${palindrome_clip_threshold}:${simple_clip_threshold}"
fi


# single-end trimming
if [ "$pe" != true ]; then
    SEQ_TYPE="SE"

    if [ -z "$path_to_trimmomatic" ] || [ "$path_to_trimmomatic" = false ]; then
                echo "Path to Trimmomatic is not set. Please set 'path_to_trimmomatic' variable in settings.txt if Trimmomatic is not installed in your PATH."
                trimmomatic $SEQ_TYPE -threads $threads -phred33 $ILLUMINACLIP \
                    "$FILE" \
                    "${main_folder}/trimmed_data/${BASE_NAME}.fq.gz" \
                    $crop $headcrop \
                    $sliding_window $leading $trailing $min_length -trimlog "${main_folder}/trimmed_data/${BASE_NAME}.trimlog"
            else
                echo "Using Trimmomatic at: $path_to_trimmomatic"
                java -jar "$path_to_trimmomatic" $SEQ_TYPE -threads $threads -phred33 $ILLUMINACLIP\
                    "$FILE" \
                    "${main_folder}/trimmed_data/${BASE_NAME}.fq.gz" \
                    $crop $headcrop \
                    $sliding_window $leading $trailing $min_length -trimlog "${main_folder}/trimmed_data/${BASE_NAME}.trimlog"
            fi
else
    SEQ_TYPE="PE"

    # Get all unique base names (e.g., name1, name2)
    # We sort and use 'uniq' to ensure each base name is processed once
    BASENAMES=$(ls -1 "${main_folder}/raw_data"/*_1.fq.gz 2>/dev/null | \
            xargs -n 1 basename | \
            sed 's/_1.fq.gz//' | sort | uniq)

    # Loop through each unique base name
    for BASE_NAME in $BASENAMES; do
        FILE1="${main_folder}/raw_data/${BASE_NAME}_1.fq.gz"
        FILE2="${main_folder}/raw_data/${BASE_NAME}_2.fq.gz"

        # Check if both paired files exist
        if [ -f "$FILE1" ] && [ -f "$FILE2" ]; then
            echo "Processing pair: $FILE1 and $FILE2"

            if [ -z "$path_to_trimmomatic" ] || [ "$path_to_trimmomatic" = false ]; then
                echo "Path to Trimmomatic is not set. Please set 'path_to_trimmomatic' variable in settings.txt if Trimmomatic is not installed in your PATH."
                trimmomatic $SEQ_TYPE -threads $threads -phred33 $ILLUMINACLIP \
                    "$FILE1" "$FILE2" \
                    "${main_folder}/trimmed_data/paired/${BASE_NAME}_1P.fq.gz" \
                    "${main_folder}/trimmed_data/unpaired/${BASE_NAME}_1U.fq.gz" \
                    "${main_folder}/trimmed_data/paired/${BASE_NAME}_2P.fq.gz" \
                    "${main_folder}/trimmed_data/unpaired/${BASE_NAME}_2U.fq.gz" \
                    $crop $headcrop \
                    $sliding_window $leading $trailing $min_length -trimlog "${main_folder}/trimmed_data/${BASE_NAME}.trimlog"
            else
                echo "Using Trimmomatic at: $path_to_trimmomatic"
                java -jar "$path_to_trimmomatic" $SEQ_TYPE -threads $threads -phred33 $ILLUMINACLIP \
                "$FILE1" "$FILE2" \
                "${main_folder}/trimmed_data/paired/${BASE_NAME}_1P.fq.gz" \
                "${main_folder}/trimmed_data/unpaired/${BASE_NAME}_1U.fq.gz" \
                "${main_folder}/trimmed_data/paired/${BASE_NAME}_2P.fq.gz" \
                "${main_folder}/trimmed_data/unpaired/${BASE_NAME}_2U.fq.gz" \
                $crop $headcrop \
                $sliding_window $leading $trailing $min_length -trimlog "${main_folder}/trimmed_data/${BASE_NAME}.trimlog"
            fi
        else
            echo "Warning: Missing paired file(s) for basename '$BASE_NAME'. Skipping."
            if [ ! -f "$FILE1" ]; then
                echo "  Missing: $FILE1"
            fi
            if [ ! -f "$FILE2" ]; then
            echo "  Missing: $FILE2"
        fi
    fi
done
fi

echo "Trimming finished."