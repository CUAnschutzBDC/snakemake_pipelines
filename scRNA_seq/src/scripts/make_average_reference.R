library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(clustifyr)

setwd("/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/src/scripts")
source("functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
  clusters <- "SCT_cluster"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
  clusters <- "RNA_cluster"
}

# Set directories
base_dir <- "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/"

save_dir <- paste0(base_dir, "results/R_analysis/")

ref_dir <- "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/files/thymus_annotated_matrix_files/"

sample <- "Hashtag_STOC"

meta_cluster <- "Anno_level_fig1"

# Read in the data
stoc_data <- readRDS(paste0(save_dir, "stoc_processed.rds"))

ref_sce <- readRDS(paste0(ref_dir, "sce_object.rda"))

ref_metadata <- colData(rec_sce)
ref_counts <- assays(rec_sce)$logcounts

ref_genes <- rowData(rec_sce)

rownames(ref_counts) <- ref_genes$X
colnames(ref_counts) <- rownames(ref_metadata)
ref_metadata <- data.frame(ref_metadata)

#seurat <- CreateSeuratObject(ref_counts)
# Set the expression assay
#seurat <- SetAssayData(seurat, "data", ref_counts)
# Add observation metadata
#seurat <- AddMetaData(seurat, ref_metadata)

#Idents(seurat) <- meta_cluster

#ref_mat <- AverageExpression(seurat)

# Find Average expression by cluster
ref_mat <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster,
  if_log = TRUE
)

write.csv(ref_mat, paste0(ref_dir, "average_cluster_expression.csv"))

