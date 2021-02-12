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

sample <- "Hashtag_STOC"

# Read in the data
stoc_data <- readRDS(paste0(save_dir, "stoc_processed.rds"))

# RNA heatmap
Idents(stoc_data) <- "SCT_cluster"
marker_genes_rna <- FindAllMarkers(stoc_data, assay = normalization_method)
write.csv(marker_genes_rna, file = paste0(save_dir,
                                          "files/rna_markes_gene_clusters.csv"))
marker_genes_adt <- FindAllMarkers(stoc_data, assay = "ADT")
write.csv(marker_genes_adt, file = paste0(save_dir,
                                          "files/adt_markes_gene_clusters.csv"))
DefaultAssay(stoc_data) <- normalization_method
top10_rna <- marker_genes_rna %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
rna_heatmap <- DoHeatmap(stoc_data, features = top10_rna$gene) + NoLegend()
DefaultAssay(stoc_data) <- "ADT"
top10_adt <- marker_genes_adt %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
adt_heatmap <- DoHeatmap(stoc_data, features = top10_adt$gene) + NoLegend()

# ADT heatmap
Idents(stoc_data) <- "ADT_cluster"
marker_genes_rna <- FindAllMarkers(stoc_data, assay = normalization_method)
write.csv(marker_genes_rna, file = paste0(save_dir,
                                          "files/rna_markes_adt_clusters.csv"))
marker_genes_adt <- FindAllMarkers(stoc_data, assay = "ADT")
write.csv(marker_genes_adt, file = paste0(save_dir,
                                          "files/adt_markes_rna_clusters.csv"))
DefaultAssay(stoc_data) <- normalization_method
top10_rna <- marker_genes_rna %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
rna_heatmap <- DoHeatmap(stoc_data, features = top10_rna$gene) + NoLegend()
DefaultAssay(stoc_data) <- "ADT"
top10_adt <- marker_genes_adt %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
adt_heatmap <- DoHeatmap(stoc_data, features = top10_adt$gene) + NoLegend()


# Combined heatmap
Idents(stoc_data) <- "combined_cluster"
marker_genes_rna <- FindAllMarkers(stoc_data, assay = normalization_method)
write.csv(marker_genes_rna, file = paste0(save_dir,
                                          "files/rna_markes_combined_clusters.csv"))
marker_genes_adt <- FindAllMarkers(stoc_data, assay = "ADT")
write.csv(marker_genes_adt, file = paste0(save_dir, "files/adt_markes_combined_clusters.csv"))
DefaultAssay(stoc_data) <- normalization_method
top10_rna <- marker_genes_rna %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
rna_heatmap <- DoHeatmap(stoc_data, features = top10_rna$gene) + NoLegend()
DefaultAssay(stoc_data) <- "ADT"
top10_adt <- marker_genes_adt %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
adt_heatmap <- DoHeatmap(stoc_data, features = top10_adt$gene) + NoLegend()
