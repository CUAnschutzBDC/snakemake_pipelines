library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "sample"

vars.to.regress <- NULL

HTO <- FALSE
ADT <- FALSE
adt_PCA <- FALSE # Set to true to run PCA on ADT data (if there are enough
                 # proteins ex: > 100). Otherwise, set to FALSE and the
                 # normalized values will be used.


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
if(ADT){
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_adt.rds"))
} else {
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_doublet.rds"))
  
}
# Remove doublets
Idents(seurat_data) <- "Doublet_finder"
seurat_data <- subset(x = seurat_data, idents = "Singlet")
if(HTO){
  Idents(seurat_data) <- "HTO_classification.global"
  seurat_data <- subset(x = seurat_data, idents = "Singlet")
}

seurat_data <- seurat_data %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# PCA --------------------------------------------------------------------------

# PCA of gene expression
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                      sample_object = seurat_data)

pdf(file.path(save_dir, "images", "RNA_pca.pdf"))
print(RNA_plots)
dev.off()

if(ADT){
  if(adt_PCA){
  # PCA of surface protein
    seurat_data <- PCA_dimRed(seurat_data, assay = "ADT")

    ADT_plots <- plot_PCA(HTO = HTO, assay = "ADT", sample_object = seurat_data)

  } else {
    # set up dsb values to use in WNN analysis 
    DefaultAssay(seurat_data) <- "ADT"
    
    # hack seurat to use normalized protein values as a dimensionality reduction object.
    VariableFeatures(seurat_data) <- rownames(seurat_data)
    
    # run true pca to initialize dr pca slot for WNN 
    seurat_data <- ScaleData(seurat_data, verbose = FALSE) %>%
      RunPCA(reduction.name = "pdsb",
             features = VariableFeatures(seurat_data),
             verbose = FALSE)
    
    # make matrix of norm values to add as dr embeddings
    pseudo <- t(GetAssayData(seurat_data, slot = "data"))
    pseudo_colnames <- paste('pseudo', 1:ncol(pseudo), sep = "_")
    colnames(pseudo) <- pseudo_colnames
    # add to object 
    seurat_data@reductions$pdsb@cell.embeddings = pseudo
    
    ADT_plots <- plotDimRed(seurat_data,
                            col_by = c("orig.ident", "percent.mt",
                                       "nFeature_RNA", "nCount_RNA",
                                       "nFeature_ADT", "nCount_ADT"),
                            plot_type = "pdsb")
  }
  
  pdf(file.path(save_dir, "images", "ADT_pca.pdf"))
  print(ADT_plots)
  dev.off()
}

# I think I won't do any batch correction, the more activated from 2 samples
# line up and the less activated from 2 samples line up. To me this indicates
# that correction is unnecessary.
#seurat_data <- RunHarmony(seurat_data, "HTO_classification",
#                          plot_convergence = TRUE, lambda = 1,
#                          theta = 14)

#harmony_plot <- plotDimRed(seurat_data, col_by = "orig.ident",
#                           plot_type = "harmony")

#cell_type_plots <- plotDimRed(seurat_data, col_by = "HTO_classification",
#                              plot_type = "harmony")

# UMAP -------------------------------------------------------------------------
RNA_pcs <- 28
ADT_pcs <- 8

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
  resolution = 0.8, assay = seurat_assay, HTO = HTO)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

plotDimRed(seurat_data, col_by = "Ins1", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "Ins2", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "Gcg", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "Sst", plot_type = "rna.umap")
plotDimRed(seurat_data, col_by = "tdTomato", plot_type = "rna.umap")

if(ADT){
  if(adt_PCA){
    adt_reduction <- "apca"
  } else{
    adt_reduction <- "pdsb"
  }
  # UMAP of surface protein
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
    resolution = 0.2, assay = "ADT", HTO = HTO, reduction = adt_reduction)

  seurat_data <- umap_data$object

  adt_plots <- umap_data$plots


  plots <- featDistPlot(seurat_data, geneset = adt_list,
                       sep_by = "ADT_cluster", combine = FALSE,
                       assay = "ADT")
  
  # WNN ------------------------------------------------------------------------
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
    seurat_data, reduction.list = list(pca_slot, adt_reduction), 
    dims.list = list(1:RNA_pcs, 1:ADT_pcs),
    modality.weight.name = c(weight_name, "ADT.weight")
  )


  seurat_data <- RunUMAP(seurat_data, nn.name = "weighted.nn",
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_data <- FindClusters(seurat_data, graph.name = "wsnn",
                                algorithm = 3, resolution = 0.6, verbose = FALSE)

  seurat_data[["combined_cluster"]] <- Idents(seurat_data)
  col_by_list <- c("combined_cluster", "orig.ident", adt_list)
  if(HTO){
    col_by_list <- c(col_by_list, "HTO_classification")
  }
  save_plot <- file.path(save_dir, "images", "combined_umap.pdf")
  plot_list <- plotDimRed(sample_object = seurat_data,
                          save_plot = save_plot,
                          col_by = col_by_list,
                          plot_type = "wnn.umap")
  VlnPlot(seurat_data, features = "RNA.weight",
          group.by = 'combined_cluster',
          sort = TRUE,
          pt.size = 0.1) +
    NoLegend()
  VlnPlot(seurat_data, features = "ADT.weight",
          group.by = 'combined_cluster',
          sort = TRUE,
          pt.size = 0.1) +
    NoLegend()
}


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
