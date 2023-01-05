library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "sample_name"

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

ADT <- FALSE
HTO <- FALSE

vars.to.regress <- NULL

# Set directories
base_dir <- here()

save_dir <- file.path(base_dir, "results", sample, "R_analysis")

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

mt_pattern <- "^mt-" # "^MT-" for human, "^mt-" for mice

# Create seurat object
seurat_object <- create_seurat_object(sample = sample,
                                      count_path = file.path(base_dir,
                                                             "results"),
                                      ADT = ADT, hashtag = HTO
                                      )

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = mt_pattern)

# Quality plots to determin cutoffs
rna_qual <- VlnPlot(seurat_object,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3)

if(ADT){
  adt_qual <- VlnPlot(seurat_object, features = c("nCount_ADT", "nFeature_ADT"))
}

pdf(file.path(save_dir, "images", "quality_plots.pdf"))
rna_qual
if(ADT){
  adt_qual
}
dev.off()

# Save before moving on
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_unfilt.rds"))

# Remove outliers
if(ADT){
  seurat_object <- subset(x = seurat_object, subset = percent.mt < 10 &
                          nFeature_RNA > 1000 & nFeature_RNA < 8500 & 
                          nCount_ADT < 10000)
} else {
  seurat_object <- subset(x = seurat_object, subset = percent.mt < 10 &
                          nFeature_RNA > 500 & nFeature_RNA < 8000)
}

# Quality plots to check cutoffs
rna_qual <- VlnPlot(seurat_object,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3)

# Normalization
# Single Cell Transform normalization
seurat_object <- SCTransform(seurat_object, vars.to.regress = vars.to.regress,
                             verbose = FALSE)

# Default normalization
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

if(HTO){
  # Demultiplex
  seurat_object <- HTODemux(seurat_object, assay = "HTO", positive.quantile = 0.90)


  # Plots
  ridge_p <- RidgePlot(seurat_object, assay = "HTO",
                       features = rownames(stoc_data[["HTO"]])[1:2], ncol = 2)
  scatter_p <- FeatureScatter(seurat_object,
                              feature1 = "hto_Hashtag-STOC-86-008",
                              feature2 = "hto_Hashtag-STOC-86-009")
  Idents(seurat_object) <- "HTO_classification.global"
  vln_p <- VlnPlot(seurat_object, features = "nCount_RNA",
                   pt.size = 0.1, log = TRUE)

  pdf(file.path(save_dir, "images", "HTO_plots.pdf"))
  ridge_p
  scatter_p
  vln_p
  dev.off()
}

saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))