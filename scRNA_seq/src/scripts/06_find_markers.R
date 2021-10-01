library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(openxlsx)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "sample"
cell_types <- "combined_celltype"
clusters <- "combined_cluster"

normalization_method <- "log" # can be SCT or log

HTO <- FALSE
ADT <- FALSE # Set to true if you want to run DE on ADT (not enough ADT here)
pval <- 0.05
logfc <- 0.5

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- "path/to/base/dir"

source(file.path(base_dir, "src", "scripts", "functions.R"))

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

gene_lists <- NULL

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Cell type DE -----------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = cell_types,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = cell_types,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# RNA cluster DE ---------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = clusters,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = clusters,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}


# Activation -------------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = "activation",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = "activation",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

seurat_data[[paste0(cell_types, sample)]] <- 
  seurat_data[[cell_types]]

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
