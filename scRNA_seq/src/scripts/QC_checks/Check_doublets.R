library(Seurat)
library(LaCroixColoR)
setwd("/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi")
source("src/scripts/functions.R")

processed_obj <- readRDS("results/R_analysis/stoc_processed.rds")

doublet_obj <- readRDS("results/R_analysis/stoc_doublet.rds")

cell_types <- processed_obj$combined_celltype

doublet_obj <- AddMetaData(doublet_obj, cell_types,
                           col.name = "combined_celltype")

# Colors
colors <- lacroix_palette("Coconut", 6, "discrete")

colors <- colors[c(1,2,3,5)]

names(colors) <- c("DC2", "Fb_cycling", "Mast", "mTEC.III.")

plotDimRed(doublet_obj, "combined_celltype", plot_type = "rna.umap",
           color = colors)

plotDimRed(doublet_obj, "DF.classifications_0.25_0.005_120",
           plot_type = "rna.umap")

plotDimRed(doublet_obj, "HTO_classification.global",
           plot_type = "rna.umap")
