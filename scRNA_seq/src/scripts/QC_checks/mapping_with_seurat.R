library(Seurat)
library(tidyverse)
library(cowplot)
library(SingleCellExperiment)
source("src/scripts/functions.R")
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

stoc_object <- readRDS("results/R_analysis/stoc_processed.rds")

thymic_reference_sce <- readRDS("files/thymus_annotated_matrix_files/sce_object.rda")

rownames(thymic_reference_sce) <- rowData(thymic_reference_sce)$X

exprs <- assays(thymic_reference_sce)$logcounts

rownames(exprs) <- rownames(thymic_reference_sce)
colnames(exprs) <- colnames(thymic_reference_sce)

metadata <- colData(thymic_reference_sce)

umap <- reducedDim(thymic_reference_sce)

# This fails with the error:
# Error in asMethod(object) : 
#Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102
# I likely need to run on the server
thymic_reference_seurat <- CreateSeuratObject(exprs)

thymic_reference_seurat <- AddMetadata(colData(thymic_reference_sce))
