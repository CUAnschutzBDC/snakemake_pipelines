library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "merge_all"

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

objects <- c("sample_1", "sample_2")

obj_list <- lapply(objects, function(x){
  
  seurat_object_path <- file.path(base_dir, "results",
                                  x, "R_analysis/rda_obj/seurat_doublet.rds")
  print(seurat_object_path)
  seurat_object <- readRDS(seurat_object_path)
  
  Idents(seurat_object) <- "Doublet_finder"
  seurat_object <- subset(x = seurat_object, idents = "Singlet")
  return(seurat_object)
})

# Merge Seurat objects
seurat_data <- merge(obj_list[[1]], obj_list[2:length(obj_list)])

seurat_data$sample <- seurat_data$orig.ident

seurat_data$individual <- gsub("_.*", "", seurat_data$sample)
seurat_data$treatment <- gsub(".*_", "", seurat_data$sample)

# Normalization
seurat_data <- seurat_data %>%
  FindVariableFeatures() %>%
  ScaleData()

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                      "seurat_start.rds"))