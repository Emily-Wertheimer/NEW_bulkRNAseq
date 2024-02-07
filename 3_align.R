## align_single_sample.R
# this is what is called by the shell scipt 

# picks up sample name from shell script loop
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

# load library for alignment
library(Rsubread)

# Define paths
fastq_dir <- "/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/sample_dir/merged_all_lanes_old"
output_dir <- "/gpfs/gibbs/pi/huckins/ekw28/ghre_lep_sema_RNAseq/bam_aligned_files"
index_path <- '/gpfs/gibbs/pi/huckins/ekw28/running_alignment_jan27/new_index/hg38_index2' 


# Define file paths
readfile1 <- file.path(fastq_dir, paste0(sample_name, "_R1_all_lanes.fastq.gz"))
readfile2 <- file.path(fastq_dir, paste0(sample_name, "_R2_all_lanes.fastq.gz"))
output_file <- file.path(output_dir, paste0(sample_name, ".BAM"))

align_and_save <- function(sample_name) {  

  # Define file paths
  readfile1 <- file.path(fastq_dir, paste0(sample_name, "_R1_all_lanes.fastq.gz"))
  readfile2 <- file.path(fastq_dir, paste0(sample_name, "_R2_all_lanes.fastq.gz"))
  output_file <- file.path(output_dir, paste0(sample_name, ".BAM"))

  # Run the alignment
  align(index = index_path,
        readfile1 = readfile1,
        readfile2 = readfile2,
        output_file = output_file,
        type = 'rna',
        minFragLength = 50,
        maxFragLength = 600
       )
}

align_and_save(sample_name)

