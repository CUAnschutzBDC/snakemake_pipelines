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

sample <- "Hashtag_STOC"

meta_cluster <- "Anno_level_fig1"

HTO <- TRUE

# Read in the data
stoc_data <- readRDS(paste0(save_dir, "stoc_processed.rds"))
DefaultAssay(stoc_data) <- seurat_assay
# Information for cell mapping
ref_dir <-
  "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/files/thymus_annotated_matrix_files/"

ref_mat <- read.csv(paste0(ref_dir, "average_cluster_expression.csv"),
                    header = TRUE, row.names = 1)

ref_mat <- ref_mat[rownames(ref_mat) %in% rownames(stoc_data),]

##########
# Set up #
##########
# get count matrix
DefaultAssay(stoc_data) <- seurat_assay
stoc_mat <- GetAssayData(object = stoc_data, slot = "data")
stoc_mat <- as.data.frame(stoc_mat)

# Find only 1000 variable features

stoc_var <- FindVariableFeatures(
  stoc_data,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 1000
)
stoc_genes <- VariableFeatures(stoc_var)

# RNA
if(SCT){
  Idents(stoc_data) <- "SCT_cluster"
  clusters = "SCT_cluster"
} else{
  Idents(stoc_data) <- "RNA_cluster"
  clusters = "RNA_cluster"
}


stoc_metadata <- stoc_data[[clusters]]

stoc_metadata[[clusters]] <- as.character(stoc_metadata[[clusters]])

stoc_data[[clusters]] <- stoc_metadata[[clusters]]

# Run clustify
stoc_res <- clustify(
  input = stoc_mat,
  metadata = stoc_metadata,
  ref_mat = ref_mat,
  query_genes = stoc_genes,
  cluster_col = clusters
)

stoc_cluster <- cor_to_call(stoc_res)
new_clusters <- stoc_cluster$type
names(new_clusters) <- stoc_cluster$cluster

stoc_data$RNA_celltype <- new_clusters[stoc_data$RNA_cluster]

plot1 <- plotDimRed(stoc_data, col_by = "RNA_celltype", plot_type = "rna.umap")
plot2 <- plotDimRed(stoc_data, col_by = "RNA_celltype", plot_type = "adt.umap")
plot3 <- plotDimRed(stoc_data, col_by = "RNA_celltype", plot_type = "wnn.umap")

pdf(paste0(save_dir, "images/RNA_celltype_mapping.pdf"))
plot1
plot2
plot3
dev.off()

# ADT
clusters <- "ADT_cluster"
stoc_metadata <- stoc_data[[clusters]]

stoc_metadata[[clusters]] <- as.character(stoc_metadata[[clusters]])

stoc_data[[clusters]] <- stoc_metadata[[clusters]]

# Run clustify
stoc_res <- clustify(
  input = stoc_mat,
  metadata = stoc_metadata,
  ref_mat = ref_mat,
  query_genes = stoc_genes,
  cluster_col = clusters
)

stoc_cluster <- cor_to_call(stoc_res)
new_clusters <- stoc_cluster$type
names(new_clusters) <- stoc_cluster$cluster

stoc_data$ADT_celltype <- new_clusters[stoc_data$ADT_cluster]

plot1 <- plotDimRed(stoc_data, col_by = "ADT_celltype", plot_type = "rna.umap")
plot2 <- plotDimRed(stoc_data, col_by = "ADT_celltype", plot_type = "adt.umap")
plot3 <- plotDimRed(stoc_data, col_by = "ADT_celltype", plot_type = "wnn.umap")

pdf(paste0(save_dir, "images/ADT_celltype_mapping.pdf"))
plot1
plot2
plot3
dev.off()

# Combined
clusters <- "combined_cluster"
stoc_metadata <- stoc_data[[clusters]]

stoc_metadata[[clusters]] <- as.character(stoc_metadata[[clusters]])

stoc_data[[clusters]] <- stoc_metadata[[clusters]]

# Run clustify
stoc_res <- clustify(
  input = stoc_mat,
  metadata = stoc_metadata,
  ref_mat = ref_mat,
  query_genes = stoc_genes,
  cluster_col = clusters
)

stoc_cluster <- cor_to_call(stoc_res)
new_clusters <- stoc_cluster$type
names(new_clusters) <- stoc_cluster$cluster

stoc_data$combined_celltype <- new_clusters[stoc_data$combined_cluster]

plot1 <- plotDimRed(stoc_data, col_by = "combined_celltype", plot_type = "rna.umap")
plot2 <- plotDimRed(stoc_data, col_by = "combined_celltype", plot_type = "adt.umap")
plot3 <- plotDimRed(stoc_data, col_by = "combined_celltype", plot_type = "wnn.umap")

pdf(paste0(save_dir, "images/combined_celltype_mapping.pdf"))
plot1
plot2
plot3
dev.off()