#!/bin/bash

# extract fastQs
base_url="http://fcb.ycga.yale.edu:3010/R81qz3tdRb85MQriRCiTqsaplcr38se/sample_dir_000015462"
sample_list_file="sample_list.txt"
OUTDIR="/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/sema_dose_curve_0924/0_raw_data"

# Loop through each sample in the list
while read -r sample; do
    # Define the sample URL
    sample_url="${base_url}/${sample}/"
    
    # recursively download files in the sample's subdirectory
    wget -r -np -nH --cut-dirs=3 -R "index.html*" -P "OUTDIR" "$sample_url"
    
done < "$sample_list_file"

# rearrange files so that all .gz fastQs are in one folder
source_dir="/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/sema_dose_curve_0924/0_raw_data/data/ycga_raw_files"
target_dir="/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/sema_dose_curve_0924/0_raw_data/data/YCGA_raw_fastQ"

find "$source_dir" -name "*.fastq.gz" -exec mv {} "$target_dir" \;


