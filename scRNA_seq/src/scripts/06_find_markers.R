library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "sample"

cell_types <- "RNA_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

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

seurat_data$cluster_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[cell_types]][[1]])

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = "cluster_celltype",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = "cluster_celltype",
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}


# Tomato+ vs - -----------------------------------------------------------------
plot(featDistPlot(seurat_data, geneset = "tdTomato", sep_by = "RNA_celltype"))

Idents(seurat_data) <- "RNA_celltype"
seurat_alpha <- subset(seurat_data, idents = "beta")

plot(featDistPlot(seurat_alpha, geneset = "tdTomato",
                  sep_by = "RNA_cluster"))


seurat_alpha$tdTomato_expr <- "negative"
seurat_alpha$tdTomato_expr[GetAssayData(seurat_alpha, slot = "data")["tdTomato",]
                           > 0.15] <- "positive"


plot(featDistPlot(seurat_alpha, geneset = "tdTomato",
                  sep_by = "tdTomato_expr"))

marker_list <- find_write_markers(seurat_object = seurat_alpha,
                                  meta_col = "tdTomato_expr",
                                  pval = pval,
                                  logfc = logfc,
                                  assay = "RNA",
                                  save_dir = save_dir)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
