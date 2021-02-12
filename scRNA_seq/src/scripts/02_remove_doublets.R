library(DoubletFinder)
library(Seurat)
save_dir <- "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/results/R_analysis/"
stoc_path <- paste0(save_dir, "stoc_start.rds")
stoc_data <- readRDS(stoc_path)

pK <- 0.04
SCT <- TRUE
seurat_assay <- "RNA"

source("functions.R")

# Remove "negatives"
Idents(stoc_data) <- "HTO_classification.global"
stoc_data <- subset(x = stoc_data, idents = c("Singlet", "Doublet"))
DefaultAssay(stoc_data) <- seurat_assay
stoc_data <- PCA_dimRed(stoc_data, assay = seurat_assay)
umap_data <- group_cells(stoc_data, assay = seurat_assay)

stoc_data <- umap_data[[1]]


## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.stoc_data <- paramSweep_v3(stoc_data, PCs = 1:10, sct = SCT)
gt.calls <- stoc_data@meta.data[rownames(sweep.res.stoc_data[[1]]),
                                    "HTO_classification.global"]
gt.calls <- as.factor(gt.calls)
sweep.stats_organoid <- summarizeSweep(sweep.res.stoc_data, GT = TRUE,
                                       GT.calls = gt.calls)
bcmvn_organoid <- find.pK(sweep.stats_organoid)

max_AUC <- max(bcmvn_organoid$MeanAUC)
max_BC <- max(bcmvn_organoid$MeanBC)
max_BCmetric <- max(bcmvn_organoid$BCmetric)
# Could use pK = 0.08 or pk = 0.19
# 0.08 is closer to the max of the AUC
pK <- 0.005
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- stoc_data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.039*nrow(stoc_data@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Find pK from bcmvn output, pN selection is less important----------------------------------------------------------------

stoc_data <- doubletFinder_v3(stoc_data, PCs = 1:10, pN = 0.25,
                                  pK = pK, nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = SCT)



saveRDS(stoc_data, paste0(save_dir, "stoc_doublet.rds"))
