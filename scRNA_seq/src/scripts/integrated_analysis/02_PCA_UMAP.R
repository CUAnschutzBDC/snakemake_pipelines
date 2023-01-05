library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(clustree)
library(harmony)
library(singlecellmethods)
library(SeuratWrappers)


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

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj/seurat_start.rds"))

# PCA --------------------------------------------------------------------------


# PCA of gene expression, weighted by cell number per sample
seurat_data <- singlecellmethods::RunBalancedPCA(obj = seurat_data, 
                                                 weight.by = "orig.ident",
                                                 npcs = 50)

RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                      sample_object = seurat_data,
                      jackstraw = FALSE)

pdf(file.path(save_dir, "images/RNA_pca.pdf"))
RNA_plots
dev.off()

if(ADT){
  # PCA of surface protein
  seurat_data <- PCA_dimRed(seurat_data, assay = "ADT")
  
  ADT_plots <- plot_PCA(HTO = HTO, assay = "ADT", sample_object = seurat_data)
  
  pdf(paste0(save_dir, "images/RNA_pca.pdf"))
  plot(ADT_plots)
  dev.off()
}

# UMAP -------------------------------------------------------------------------

RNA_pcs <- 30
ADT_pcs <- 8

seurat_data$RNA_cluster <- NULL

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data <- FindClusters(seurat_data, resolution = c(0.2, 0.5, 0.8, 1, 1.2))
clustree(seurat_data)

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.8, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

plotDimRed(seurat_data, col_by = "RNA_celltype", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "sample", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "group", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "batch", plot_type = "rna.umap")

seurat_data <- BuildClusterTree(seurat_data, dims = 1:RNA_pcs)
PlotClusterTree(seurat_data)

seurat_data$uncorrected_cluster <- seurat_data$RNA_cluster

seurat_data$RNA_cluster <- NULL

# Batch correction -------------------------------------------------------------

# NOTE: Decide if batch correction is needed
# NOTE: Test Harmony and MNN to decide what is needed.

## Harmony ---------------------------------------------------------------------

seurat_data <- RunHarmony(seurat_data, c("batch"),
                          plot_convergence = TRUE)

harmony_plot <- plotDimRed(seurat_data, col_by = "sample",
                           plot_type = "harmony")

## MNN -------------------------------------------------------------------------
sce_data <- as.SingleCellExperiment(seurat_data)

set.seed(0)
# Run fastMNN on batch
corrected_data <- fastMNN(sce_data, batch = sce_data$orig.ident)

# Check reduced dims
if(!identical(rownames(SingleCellExperiment::reducedDim(x = corrected_data)), 
              colnames(seurat_data))){
  stop("names are not the same between the object and fastmnn output!!!")
}

# Change names
colnames(SingleCellExperiment::reducedDim(x = corrected_data)) <-
  paste0("mnn_", 1:ncol(SingleCellExperiment::reducedDim(x = corrected_data)))

# Add to seurat object
seurat_data$mnn <- CreateDimReducObject(
  embeddings = SingleCellExperiment::reducedDim(x = corrected_data),
  loadings = as.matrix(SingleCellExperiment::rowData(x = corrected_data)),
  assay = DefaultAssay(object = seurat_data),
  key = "mnn_"
)

# UMAP -------------------------------------------------------------------------

reduction_use <- "mnn" # can be "harmony" or "mnn"
RNA_pcs <- 20
ADT_pcs <- 8

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}

## Harmony ---------------------------------------------------------------------
set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = HTO,
                         reduction = "harmony")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data$harmony_clust <- seurat_data$RNA_cluster

cm <- confusionMatrix(seurat_data$harmony_clust, seurat_data$RNA_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

## MNN -------------------------------------------------------------------------

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = npcs,
                         resolution = 0.6, assay = seurat_assay, HTO = HTO,
                         reduction = "mnn", size = 1)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

seurat_data$mnn_clust <- seurat_data$RNA_cluster

cm <- confusionMatrix(seurat_data$mnn_clust, seurat_data$RNA_celltype)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

## Compare --------------------------------------------------------------------

cm <- confusionMatrix(seurat_data$harmony_clust, seurat_data$mnn_clust)

cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

plot1 <- plotDimRed(seurat_data, col_by = c("mnn_clust", "harmony_clust"), 
                    plot_type = "harmony.umap")

plot2 <- plotDimRed(seurat_data, col_by = c("mnn_clust", "harmony_clust"), 
                    plot_type = "mnn.umap")

## Cluster --------------------------------------------------------------------

# Pick mnn or harmony
reduction_use <- "mnn"

seurat_data <- FindClusters(seurat_data, resolution = c(0.5, 0.8, 1, 1.2))
clustree(seurat_data)

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 1, assay = seurat_assay, HTO = HTO,
                         reduction = reduction_use)

table(paste0(seurat_data$RNA_cluster, "_", seurat_data$RNA_celltype))

plot_type <- paste0(reduction_use, ".umap")

plotDimRed(seurat_data, col_by = "RNA_celltype", plot_type = plot_type)
plotDimRed(seurat_data, col_by = "sample", plot_type = plot_type)
plotDimRed(seurat_data, col_by = "group", plot_type = plot_type)
plotDimRed(seurat_data, col_by = "batch", plot_type = plot_type)

seurat_data <- BuildClusterTree(seurat_data, dims = 1:RNA_pcs)
PlotClusterTree(seurat_data)

seurat_data$uncorrected_cluster <- seurat_data$RNA_cluster

if(ADT){
  # UMAP of surface protein
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
                           resolution = 0.6, assay = "ADT", HTO = TRUE)
  
  seurat_data <- umap_data$object
  
  adt_plots <- umap_data$plots
  
  
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
  seurat_data <- FindMultiModalNeighbors(
    seurat_data, reduction.list = list(pca_slot, "apca"), 
    dims.list = list(1:RNA_pcs, 1:ADT_pcs),
    modality.weight.name = c(weight_name, "ADT.weight")
  )
  
  
  seurat_data <- RunUMAP(seurat_data, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_data <- FindClusters(seurat_data, graph.name = "wsnn",
                              algorithm = 3, resolution = 2, verbose = FALSE)
  
  seurat_data[["combined_cluster"]] <- Idents(seurat_data)
  col_by_list <- c("combined_cluster", "orig.ident")
  if(HTO){
    col_by_list <- c(col_by_list, "HTO_classification")
  }
  save_plot <- file.path(save_dir, "images/combined_umap.pdf")
  plot_list <- plotDimRed(sample_object = seurat_data,
                          save_plot = NULL,
                          col_by = col_by_list, return_plot = TRUE,
                          plot_type = "wnn.umap")
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))
