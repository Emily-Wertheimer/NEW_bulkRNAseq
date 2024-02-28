### SVA detects batch effectS and unknown noise (method 2)

# --------------------------------------------------------------------------- IMPORT LIBRARIES
library(sva)
library(limma)
library(pheatmap)

# --------------------------------------------------------------------------- IMPORT  
setwd('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG')

# metas
meta.iGABA <- read_csv('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/iGABA_meta_for_SVA.csv')
meta.iGlut <- read_csv('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/iGlut_meta_for_SVA.csv')
meta.iAstro <- read_csv('/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/iAstro_meta_for_SVA.csv')

# DGE objects

# gene counts after cpm normalization

# --------------------------------------------------------------------------- ADD BATCH AND NORM FACTORS
# add batch, lib size, norm factors cols, etc. as needed

meta.iGABA$lib.size <- data.dge.iGABA$samples$lib.size/1000000
meta.iGlut$lib.size <- data.dge.iGlut$samples$lib.size/1000000
meta.iAstro$lib.size <- data.dge.iAstro$samples$lib.size/1000000

meta.iGABA$norm.factors <- data.dge.iGABA$samples$norm.factors
meta.iGlut$norm.factors <- data.dge.iGlut$samples$norm.factors
meta.iAstro$norm.factors <- data.dge.iAstro$samples$norm.factors

# --------------------------------------------------------------------------- RUN SVA F(X)
runSVA <- function(data.dge, meta.df, gene_counts_aftercpm.df, output_prefix) {
  
  # Create models
  mod <- model.matrix(~ Treatment + lib.size , meta.df)
  mod0 <- model.matrix(~ lib.size, meta.df)
  
  # Run sva to determine the number of surrogate variables
  n.sv <- num.sv(gene_counts_aftercpm.df, mod, method="leek", B = 20)
  print(paste("Estimated number of SVs:", n.sv))
  
  # Estimate surrogate variables
  svobj <- svaseq(as.matrix(gene_counts_aftercpm.df), mod, mod0, n.sv=n.sv)
  
  # Include SVs into metadata
  svobj$sv <- as.data.frame(svobj$sv)
  meta_sv.df <- cbind(meta.df, svobj$sv)
  
  # Save the enhanced metadata with SVs
  output_path <- paste0("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/", output_prefix, "_meta_with_SVs.csv")
  write.csv(meta_sv.df, output_path, row.names = FALSE)
  
  ## Reset formula with SVs included
  form <- as.formula(paste("~ lib.size", paste(colnames(svobj$sv), collapse=" + "), sep=" + "))
  
  # Plot correlation matrix with surrogate variables
  C <- canCorPairs(form, meta_sv.df)
  graphics.off()
  p <- pheatmap(
    C, 
    color = hcl.colors(50, "YlOrRd", rev = TRUE),
    fontsize = 8,
    border_color = "black",
    cellwidth = unit(0.4, "cm"),
    cellheight = unit(0.4, "cm")
  )
  
  # Save the plot
  plot_path <- paste0("/gpfs/gibbs/pi/huckins/ekw28/bulkRNAseq/ghre_lep_sema_RNAseq/3_DEG/", output_prefix, "_correlation_matrix.png")
  png(plot_path)
  print(p)
  dev.off()
}



# --------------------------------------------------------------------------- RUN SVA

runSVA(data.dge.iAstro, meta.iAstro, gene_counts_aftercpm.df.iAstro, 'iAstro_SVA')
runSVA(data.dge.iGlut, meta.iGlut, gene_counts_aftercpm.df.iGlut, 'iGlut_SVA')
runSVA(data.dge.iGABA, meta.iGABA, gene_counts_aftercpm.df.iGABA, 'iGABA_SVA')
