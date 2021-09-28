library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)

source("functions.R")

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

# Set directories
base_dir <- "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/"

save_dir <- paste0(base_dir, "results/R_analysis/")

sample <- "sample"

# Create seurat object
stoc_data <- create_seurat_object(sample = sample,
                                  count_path = paste0(base_dir, "results"),
                                  ADT = TRUE, hashtag = TRUE
                                  )

# Add mitochondrial percent
stoc_data[["percent.mt"]] <- PercentageFeatureSet(stoc_data,
                                                  pattern = "^MT-")

# Quality plots to determin cutoffs
rna_qual <- VlnPlot(stoc_data,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3)

adt_qual <- VlnPlot(stoc_data, features = c("nCount_ADT", "nFeature_ADT"))

pdf(paste0(save_dir, "images/quality_plots.pdf"))
rna_qual
adt_qual
dev.off()
# Remove outliers
stoc_data <- subset(x = stoc_data, subset = percent.mt < 10 &
      nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_ADT < 10000)

# Normalization
if(SCT){
  # Single Cell Transform normalization
  stoc_data <- SCTransform(stoc_data, vars.to.regress = "percent.mt",
                           verbose = FALSE)
} else {
  # Default normalization
  DefaultAssay(stoc_data) <- "RNA"
  stoc_data <- NormalizeData(stoc_data) %>% 
    FindVariableFeatures() %>%
    ScaleData()
}

# Demultiplex
stoc_data <- HTODemux(stoc_data, assay = "HTO", positive.quantile = 0.90)


# Plots
ridge_p <- RidgePlot(stoc_data, assay = "HTO",
                     features = rownames(stoc_data[["HTO"]])[1:2], ncol = 2)
scatter_p <- FeatureScatter(stoc_data,
                            feature1 = "hto_Hashtag-STOC-86-008",
                            feature2 = "hto_Hashtag-STOC-86-009")
Idents(stoc_data) <- "HTO_classification.global"
vln_p <- VlnPlot(stoc_data, features = "nCount_RNA",
                 pt.size = 0.1, log = TRUE)

pdf(paste0(save_dir, "images/HTO_plots.pdf"))
ridge_p
scatter_p
vln_p
dev.off()

saveRDS(stoc_data, file = paste0(save_dir, "stoc_start.rds"))

