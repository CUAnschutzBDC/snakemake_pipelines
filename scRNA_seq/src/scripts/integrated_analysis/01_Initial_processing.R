library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "sample_all"

normalization_method <- "log" # can be SCT or log

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

# Make output directories
ifelse(!dir.exists(file.path(base_dir_proj)),
       dir.create(file.path(base_dir_proj)), FALSE)

ifelse(!dir.exists(file.path(save_dir)),
       dir.create(file.path(save_dir)), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

# List of objects and any batch / treatment info
objects <- list("sample1" = list("group" = "control", "batch" = "batch1"),
                "sample2" = list("group" = "mutant", "batch" = "batch1"),
                "sample3" = list("group" = "control", "batch" = "batch2"),
                "sample4" = list("group" = "mutant", "batch" = "batch2"))

obj_list <- lapply(names(objects), function(x){
  seurat_object_path <- file.path(base_dir, "results",
                                  x, "R_analysis/rda_obj/seurat_processed.rds")
  seurat_object <- readRDS(seurat_object_path)
  seurat_object$batch <- objects[[x]]$batch
  seurat_object$group <- objects[[x]]$group
  return(seurat_object)
})

# Merge Seurat objects
seurat_data <- merge(obj_list[[1]], obj_list[2:length(obj_list)])

seurat_data$sample <- seurat_data$orig.ident

seurat_data <- seurat_data %>%
  FindVariableFeatures() %>%
  ScaleData()

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))