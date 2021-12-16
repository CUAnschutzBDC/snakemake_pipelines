library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

ADT <- TRUE
HTO <- TRUE

vars.to.regress <- "percent.mt"

# Set directories
base_dir <- "path/to/base/dir"

source(file.path(base_dir, "src", "scripts", "functions.R"))

save_dir <- file.path(base_dir, "results", "R_analysis")

sample <- "sample"

mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice

# Create seurat object
seurat_object <- create_seurat_object(sample = sample,
                                      count_path = paste0(base_dir, "results"),
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
# Remove outliers
if(ADT){
  seurat_object <- subset(x = seurat_object, subset = percent.mt < 10 &
                          nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                          nCount_ADT < 10000)
} else {
  seurat_object <- subset(x = seurat_object, subset = percent.mt < 10 &
                          nFeature_RNA > 200 & nFeature_RNA < 6000)
}

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

saveRDS(seurat_object, file = file.path(save_dir, "seurat_start.rds"))