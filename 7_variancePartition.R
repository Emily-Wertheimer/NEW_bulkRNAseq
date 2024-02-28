# --------------------------------------------------------------------------- IMPORT LIBRARIES

# install packages & dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(limma)
library(readr)
library(dplyr)
install.packages("tibble")
library(tibble)
library(edgeR)
library(data.table)
BiocManager::install("variancePartition")
### library
# basic operation
library(tidyverse)
# for plot
install.packages('ggfortify')
library(ggfortify) # PCA
library(ggplot2)
# for variance partition
library(variancePartition)
library(pheatmap)
# for SVA
BiocManager::install("sva")
library(sva)
# for RUV
BiocManager::install("RUVSeq")
library(RUVSeq)
library(RColorBrewer)

# --------------------------------------------------------------------------- VARIANCE PARTITION
# this analysis helps you select covariates
# Check covariate corr w/ PCs 
# in how many genes is expression associated w/ a given covariate (include factors (metadata.df cols) that >1% var in 10% genes)
# make sure to include in model these covariates (>1% in 10% genes); if less than that threshold then can justify removing from subsequent model

library('variancePartition')

# select only cols of interest from metadata.df
meta.iGABA <- subset(metadata.df.iGABA, select=c('sampleName', 'Donor', 'Treatment', 'Well' ))
meta.iGlut <- subset(metadata.df.iGlut, select=c('sampleName', 'Donor', 'Treatment', 'Well' ))
meta.iAstro <- subset(metadata.df.iAstro, select=c('sampleName', 'Donor', 'Treatment', 'Well' ))

# change voom mat sample names to match meta sample names
colnames(data.iGABA.voom) <- gsub("^X", "", sub("_.*", "", colnames(data.iGABA.voom)))
colnames(data.iGlut.voom) <- gsub("^X", "", sub("_.*", "", colnames(data.iGlut.voom)))
colnames(data.iAstro.voom) <- gsub("^X", "", sub("_.*", "", colnames(data.iAstro.voom)))

# model (all vars are categorical --> modeled as random effects)
form <- ~ (1|Donor) + (1|Treatment) + (1|Well) 
form.iAstro <- ~ (1|Treatment) + (1|Well) 


# 1.  fit model and extract results
# since all vars are categorical, linear model is used 
# e/ entry  in results is a regression model fit on single gene

# 2. exact variance fractions from e/ model fit 
# for e/ gene, returns fraction of variance attributable to each variable

# INTERPRETATION
#The variance explained by e/ variable after correcting for all other variables
varPart.iGABA <- fitExtractVarPartModel(data.iGABA.voom, form, meta.iGABA)
varPart.iGlut <- fitExtractVarPartModel(data.iGlut.voom, form, meta.iGlut)
varPart.iAstro <- fitExtractVarPartModel(data.iAstro.voom, form.iAstro, meta.iAstro)

# sort vars by median fraction of variance explained
vp.iGABA <- sortCols(varPart.iGABA)
vp.iGlut <- sortCols(varPart.iGlut)
vp.iAstro <- sortCols(varPart.iAstro)

# write file
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/')
write.table(vp.iGABA, sprintf("variancePartition.iGABA_output.txt"), sep="\t", quote=F, row.names=F)
write.table(vp.iGlut, sprintf("variancePartition.iGlut_output.txt"), sep="\t", quote=F, row.names=F)
write.table(vp.iAstro, sprintf("variancePartition.iAstro_output.txt"), sep="\t", quote=F, row.names=F)

################# visualize variance partition ################
library(ggplot2)
library(reshape2)

library(ggplot2)
library(reshape2)

# Function to plot variance partition for multiple sample types
plotVariancePartition <- function(vp, outputFilename = "variancePartition.png") {
  col <- c(ggColorHue(ncol(vp) - 1), "grey85")
  vp2 <- (vp * 100) # Convert to percentage
  names <- colnames(vp2)
  t <- apply(vp2, 2, FUN=median)
  t <- round(t, 2)
  t[t == 0.00] <- "<0.1"
  textMatrix <- paste(names, "\n(", t, "%)\n", sep = "")
  melted <- melt(vp2)
  
  p <- ggplot(melted, aes(x=variable, y=value, color=variable)) +
    geom_violin(aes(x = variable, y = value), scale = "width") +
    geom_boxplot(aes(x = variable, y = value), width = .1, outlier.size = 0.4) +
    theme_bw() + ylim(0, 100) +
    theme(axis.line.x = element_line(linewidth = 0.3, linetype = "solid", colour = "black"),
          axis.line.y = element_line(linewidth = 0.3, linetype = "solid", colour = "black"),
          axis.text.x = element_text(size=8, angle=45, colour="black", hjust=1),
          axis.text.y = element_text(size=7, colour="black")) +
    scale_x_discrete(labels=textMatrix) + 
    ylab("Variance Explained (%)")
  
  # Save the plot to file
  setwd("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/figures")
  ggsave(filename = outputFilename, plot = p, device = "png", width = 10, height = 7, units = "in", dpi = 300)
  
}

plotVariancePartition(vp.iGABA, "variancePartition_iGABA.png")
plotVariancePartition(vp.iGlut, "variancePartition_iGlut.png")
plotVariancePartition(vp.iAstro, "variancePartition_iAstro.png")


################# choose covariates ################
#alanna's criteria: continuous technical factors that explained ≥ 1% variation in ≥ 10% of genes
# use this info when returning to 4_DEG.R script to design models/contrasts to find DEGs

chooseCovariates <- function(vp, outputFilename = "chooseCov.txt") {
  
  cat("Choosing Covariates\n")
  test=colnames(vp)[1:length(colnames(vp))-1]
  covs_to_adjust=c()
  for(col_i in test) {
    print(col_i)
    if(length(which(vp[,col_i] >= 0.01))/(nrow(vp)) > 0.1) covs_to_adjust=c(covs_to_adjust,col_i)
  }
  
  setwd("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG")
  write.table(covs_to_adjust, file=sprintf("Covs.txt", quote=F, sep="\n", col.names=T, row.names=F))
  
}
chooseCovariates(vp.iGABA, "chooseCovariates_iGABA.txt")
chooseCovariates(vp.iGlut, "chooseCovariates_iGlut.txt")
chooseCovariates(vp.iAstro, "chooseCovariates_iAstro.txt")
