library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log
reduction_use <- "mnn" # Pick mnn or harmony
sample <- "sample_all"

all_ref_dir <-
  "/Users/wellskr/Documents/Analysis/references/single_cell_references"

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

plot_type <- paste0(reduction_use, ".umap")

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results_tomato", sample)
save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

#-------------------------------------------------------------------------------

# NOTE update the section below to include any references pertanent to your
# sample.

########################
# Pancreas development #
########################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/byrnes_2018_mouse")

ref_mat <- read.csv(file.path(ref_dir, "E14_epithelial_average.csv"),
                    header = TRUE, row.names = 1)


cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "combined_celltype_byrnes", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 1000, clusters = "RNA_cluster",
                             plot_type = plot_type)

seurat_data <- cluster_res$object

seurat_res <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_combined_celltype_byrnes",
           plot_type = plot_type)

write.csv(seurat_res, file.path(save_dir,
                                "files/celltype_mapping_byrnes_2018.csv"))


#-------------------------------------------------------------------------------

##################
# Pancreas atlas #
##################

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pancreas/Baron_2016")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_mouse_reference.csv"),
                    header = TRUE, row.names = 1)

# This looks better than trying to use the corrected cluster
cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "baron_celltype", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 1500, clusters = "RNA_cluster",
                             plot_type = plot_type)

seurat_data <- cluster_res$object

seurat_res <- cluster_res$RNA

plotDimRed(seurat_data, col_by = "RNA_baron_celltype", plot_type = plot_type)

write.csv(seurat_res,
          file.path(save_dir, "files/celltype_mapping_baron_2016.csv"))


#-------------------------------------------------------------------------------

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
