# RNA-seq_pipeline
A standard pipeline for RNA-seq analysis

requisites:
-fastQC
-multiQC (optional)
-trimmomatic
-STAR
-subread
-samtools

General usage:
1) First, fill the settings.txt file with your parameters.

2) Run setup_rnaseq_directory.sh and place your fq.gz files in raw_data folder.

3) Run main.sh to automatically compute count matrices from raw data. 
   Alternatively, you can run one step at a time to check the outputs.

