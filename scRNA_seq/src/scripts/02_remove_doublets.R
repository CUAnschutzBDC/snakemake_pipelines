library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "sample"

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

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_start.rds"))

# Remove "negatives"
if(HTO){
  Idents(seurat_data) <- "HTO_classification.global"
  seurat_data <- subset(x = seurat_data, idents = c("Singlet", "Doublet"))
}

# PCA and UMAP
set.seed(0)
DefaultAssay(seurat_data) <- seurat_assay
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)
npcs <- 20

umap_data <- group_cells(seurat_data, assay = seurat_assay, nPCs = npcs)

seurat_data <- umap_data[[1]]

# Run Doublet finder -----------------------------------------------------------
set.seed(0)
sweep.res.seurat_data <- paramSweep_v3(seurat_data, PCs = 1:npcs, sct = SCT)

if(HTO){
  ## pK Identification (ground-truth)
  gt.calls <- seurat_data@meta.data[rownames(sweep.res.seurat_data[[1]]),
                                    "HTO_classification.global"]
  gt.calls <- as.factor(gt.calls)
  sweep.stats_sample <- summarizeSweep(sweep.res.seurat_data, GT = TRUE,
                                       GT.calls = gt.calls)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  max_AUC <- max(bcmvn_sample$MeanAUC)
  max_BC <- max(bcmvn_sample$MeanBC)
  max_BCmetric <- max(bcmvn_sample$BCmetric)
} else {
  ## pK Identification (no ground-truth)
  sweep.stats_sample <- summarizeSweep(sweep.res.seurat_data, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  max_BC <- max(bcmvn_sample$MeanBC)
  max_BCmetric <- max(bcmvn_sample$BCmetric)
}

# Best is PK 0.005 second is 0.08
pK <- 0.005
## Homotypic Doublet Proportion Estimate ---------------------------------------
annotations <- seurat_data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
## CHANGE BASED ON CELL NUMBERS
# Multiplet rate	# cells loaded	# cells recovered
# 0.40%	825	500
# 0.80%	1,650	1,000
# 1.60%	3,300	2,000
# 2.40%	4,950	3,000
# 3.20%	6,600	4,000
# 4.00%	8,250	5,000
# 4.80%	9,900	6,000
# 5.60%	11,550	7,000
# 6.40%	13,200	8,000
# 7.20%	14,850	9,000
# 8.00%	16,500	10,000
nExp_poi <- round(0.032*nrow(seurat_data@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Find pK from bcmvn output, pN selection is less important--------------------

seurat_data <- doubletFinder_v3(seurat_data, PCs = 1:10, pN = 0.25,
                                  pK = pK, nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = SCT)

seurat_data$Doublet_finder <- seurat_data$DF.classifications_0.25_0.005_80

plotDimRed(seurat_data, col_by = "Doublet_finder", plot_type = "rna.umap")

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_doublet.rds"))
