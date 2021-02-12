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

HTO <- TRUE
# Set directories
base_dir <- "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/"

save_dir <- paste0(base_dir, "results/R_analysis/")

sample <- "Hashtag_STOC"


# Read in the data
stoc_data <- readRDS(paste0(save_dir, "stoc_doublet.rds"))

# Remove doublets
Idents(stoc_data) <- "DF.classifications_0.25_0.005_120"
stoc_data <- subset(x = stoc_data, idents = "Singlet")
Idents(stoc_data) <- "HTO_classification.global"
stoc_data <- subset(x = stoc_data, idents = "Singlet")

#######
# PCA #
#######

# PCA of gene expression
stoc_data <- PCA_dimRed(stoc_data, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = TRUE, assay = seurat_assay,
                      sample_object = stoc_data)
pdf(paste0(save_dir, "images/RNA_pca.pdf"))
RNA_plots
dev.off()

# PCA of surface protein
stoc_data <- PCA_dimRed(stoc_data, assay = "ADT")

ADT_plots <- plot_PCA(HTO = TRUE, assay = "ADT", sample_object = stoc_data)

pdf(paste0(save_dir, "images/RNA_pca.pdf"))
ADT_plots
dev.off()

########
# UMAP #
######## 

# UMAP of gene expression
umap_data <- group_cells(stoc_data, "STOC1", save_dir, nPCs = 27,
  resolution = 0.6, assay = seurat_assay, HTO = TRUE)

stoc_data <- umap_data$object

stoc_plots <- umap_data$plots


# UMAP of surface protein
umap_data <- group_cells(stoc_data, "STOC1", save_dir, nPCs = 16,
  resolution = 0.6, assay = "ADT", HTO = TRUE)

stoc_data <- umap_data$object

stoc_plots <- umap_data$plots


# UMAP of combined
# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
if(SCT){
  pca_slot <- "sctpca"
  weight_name <- "SCT.weight"
} else{
  pca_slot <- "pca"
  weight_name <- "RNA.weight"
}
stoc_data <- FindMultiModalNeighbors(
  stoc_data, reduction.list = list(pca_slot, "apca"), 
  dims.list = list(1:27, 1:16),
  modality.weight.name = c(weight_name, "ADT.weight")
)


stoc_data <- RunUMAP(stoc_data, nn.name = "weighted.nn",
              reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
stoc_data <- FindClusters(stoc_data, graph.name = "wsnn",
                              algorithm = 3, resolution = 2, verbose = FALSE)

stoc_data[["combined_cluster"]] <- Idents(stoc_data)
col_by_list <- c("combined_cluster", "orig.ident")
if(HTO){
  col_by_list <- c(col_by_list, "HTO_classification")
}
save_plot <- paste0(save_dir, "images/combined_umap.pdf")
plot_list <- plotDimRed(sample_object = stoc_data,
                        save_plot = NULL,
                        col_by = col_by_list, return_plot = TRUE,
                        plot_type = "wnn.umap")

saveRDS(stoc_data, paste0(save_dir, "stoc_processed.rds"))
