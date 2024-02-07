## map GeneID to actual gene name (using R)

# import libraries
library(tidyverse)
library(data.table)

# load raw counts data 
counts <- data.frame(fread("RawCountMatrix_final.txt"))

######## option 1: get gene names using Carina's file ########
anno <- read.csv("/gpfs/gibbs/pi/huckins/ekw28/resources")  # message Carina for this file if you don't already have it

######## option 2: get gene names using biomaRt ########
	library(biomaRt)
	ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") #mirror="useast" if needed

	my_genes <- meta1$MarkerName
	my_biomart <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description','chromosome_name','start_position','end_position'),
		      filters = 'ensembl_gene_id',
		      values = my_genes, 
		      mart = ensembl)

######## once you have gene names: ########

anno <- anno[,c(10,14)]
colnames(anno) <- c("Gene_name", "GeneID")
counts <- merge(anno, counts, by="GeneID")
write.table(counts, "Count_Matrix_Annotated.txt", quote=F, row.names=F)
