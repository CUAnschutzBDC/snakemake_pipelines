library(DoubletFinder)
library(Seurat)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "sample"

normalization_method <- "log" # can be SCT or log

HTO <- TRUE
ADT <- TRUE

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
umap_data <- group_cells(seurat_data, assay = seurat_assay)

seurat_data <- umap_data[[1]]

# Run Doublet finder -----------------------------------------------------------
set.seed(0)
sweep.res.seurat_data <- paramSweep_v3(seurat_data, PCs = 1:10, sct = SCT)

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

# Best is PK 0.15
pK <- 0.005
## Homotypic Doublet Proportion Estimate ---------------------------------------
annotations <- seurat_data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
## CHANGE BASED ON CELL NUMBERS
nExp_poi <- round(0.076*nrow(seurat_data@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Find pK from bcmvn output, pN selection is less important--------------------

seurat_data <- doubletFinder_v3(seurat_data, PCs = 1:10, pN = 0.25,
                                  pK = pK, nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = SCT)

seurat_data$Doublet_finder <- seurat_data$DF.classifications_0.25_0.005_700

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_doublet.rds"))
