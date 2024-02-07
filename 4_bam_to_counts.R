## takes .BAM files from alignment and makes count files for each sample

# load required libraries
library(Rsubread)
library(dplyr)

# set the R script to take arguments from the command line
args <- commandArgs(trailingOnly=TRUE)
sample_id <- args[1] # Directly use the argument as the sample identifier

project_dir <- "/gpfs/gibbs/pi/huckins/ekw28/running_alignment_jan27/bam_files/"
output_dir <- "/gpfs/gibbs/pi/huckins/ekw28/counts/count_matrices/"

bam_file_name <- sample_id #paste0(sample_id, ".BAM") # Construct the BAM file name using the sample ID

# Check if the BAM file actually exists to avoid running into the same error
if (!file.exists(paste0(project_dir, bam_file_name))) {
  stop("BAM file does not exist: ", paste0(project_dir, bam_file_name))
}

# Check if the count matrix for this sample has already been generated
output_file_name <- paste0(sample_id, "_count_matrix.txt")
if (file.exists(paste0(output_dir, output_file_name))) {
  cat("Count matrix already exists for sample: ", sample_id, "\n")
} else {
  # Perform featureCounts
  count_matrix <- featureCounts(files = paste0(project_dir, bam_file_name),
                                annot.inbuilt = "hg38",
                                isPairedEnd = TRUE,
                                useMetaFeatures = TRUE,
                                allowMultiOverlap = FALSE,
                                largestOverlap = TRUE,
                                countMultiMappingReads = TRUE)

  # Extract counts from output
  counts_df <- data.frame(GeneID = count_matrix$annotation$GeneID,
                          Length = count_matrix$annotation$Length,
                          Counts = count_matrix$counts)

  # Save the counts
  write.table(counts_df, 
              file = paste0(output_dir, output_file_name), 
              row.names = FALSE, 
              quote = FALSE, 
              sep = "\t")
  cat("Count matrix generated for sample: ", sample_id, "\n")
}
