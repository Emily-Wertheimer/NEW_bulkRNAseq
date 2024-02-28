### SVA detects batch effectS and unknown noise (method 2)

# --------------------------------------------------------------------------- IMPORT LIBRARIES

library(sva)
library(limma)
library(pheatmap)

# --------------------------------------------------------------------------- MAKE META FOR E/ CELL TYPE

# meta.iGABA <- subset(metadata.df.iGABA, select=c('Donor', 'Treatment', 'Well' ))
# meta.iGlut <- subset(metadata.df.iGlut, select=c('Donor', 'Treatment', 'Well' ))
# meta.iAstro <- subset(metadata.df.iAstro, select=c('Treatment', 'Well' )) # only have 1 donor
# add batch col
#meta.iGABA$batch <- 1
#meta.iGlut$batch <- 1
#meta.iAstro$batch <- 2

## use combined meta, add batch, lib size, norm factors cols
metadata.df$batch <- NA
metadata.df$batch[metadata.df$cellType == 'iGABA' | metadata.df$cellType == 'iGlut'] <- 1
metadata.df$batch[metadata.df$cellType == 'iAstro'] <- 2


#### should i be running this in each cell type separately? 
metadata.df$lib.size <- data.dge.all$samples$lib.size/1000000
metadata.df$norm.factors <- data.dge.all$samples$norm.factors

# --------------------------------------------------------------------------- RUN SVA F(X)

runSVA <- function(data.dge, meta.df, gene_counts_aftercpm.df) {
  
  # Update metadata with library size and normalization factors
  #meta.df$lib.size <- data.dge[["samples"]][["lib.size"]]/1000000
  #meta.df$norm.factors <- data.dge[["samples"]][["norm.factors"]]
  

  # Create models
  mod <- model.matrix(~ Treatment + batch + lib.size + norm.factors + Well, meta.df)
  mod0 <- model.matrix(~ batch + lib.size + norm.factors, meta.df)
  
  # Run sva to determine the number of surrogate variables
  n.sv <- num.sv(gene_counts_aftercpm.df, mod, method="leek", B = 20)
  print(head(n.sv))
  
  # Estimate surrogate variables
  svobj <- svaseq(as.matrix(gene_counts_aftercpm.df), mod, mod0, n.sv=4) # Adjust n.sv based on previous results if necessary
  
  # Include SVs into metadata
  svobj$sv <- as.data.frame(svobj$sv)
  meta_sv.df <- cbind(meta.df, svobj$sv)
  
  # Reset formula with SVs included
  form <- as.formula(paste("~ group + batch + lib.size + norm.factors", paste(colnames(svobj$sv), collapse=" + "), sep=" + "))
  
  # Plot correlation matrix with surrogate variables
  C <- canCorPairs(form, meta_sv.df)
  graphics.off()
  pheatmap(
    C, 
    color = hcl.colors(50, "YlOrRd", rev = TRUE),
    fontsize = 8,
    border_color = "black",
    cellwidth = unit(0.4, "cm"),
    cellheight = unit(0.4, "cm")
  )
}

# --------------------------------------------------------------------------- RUN SVA

runSVA(data.dge.iAstro, meta.iAstro, gene_counts_aftercpm.df.iAstro)
runSVA(data.dge.iGlut, meta.iGlut, gene_counts_aftercpm.df.iGlut)
runSVA(data.dge.iGABA, meta.iGABA, gene_counts_aftercpm.df.iGABA)
