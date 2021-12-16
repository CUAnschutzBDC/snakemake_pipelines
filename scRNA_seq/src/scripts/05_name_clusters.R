library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(clustifyr)
library(viridis)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "sample"

normalization_method <- "log" # can be SCT or log

HTO <- TRUE
ADT <- TRUE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
  clusters <- "SCT_cluster"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
  clusters <- "RNA_cluster"
}

# Set directories
base_dir <- "path/to/base/dir"

source(file.path(base_dir, "src", "scripts", "functions.R"))

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

## Single cell reference -------------------------------------------------------
ref_dir <- file.path(base_dir, "files", "celltype_refs", "hca_pbmc")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
DefaultAssay(seurat_data) <- seurat_assay

# Ref info
ref_obj <- readRDS(file.path(ref_dir, "processed_obj.rda"))

# Make average reference for celltype
Idents(ref_obj) <- "Celltype_Annotation"
ref_mat <- AverageExpression(ref_obj, "gene_RNA")
ref_mat <- ref_mat$gene_RNA

seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat,
                                 save_name = "celltype",
                                 ADT = ADT,
                                 nfeatures = 2000,
                                 clusters = clusters)


seurat_data <- seurat_res_list$object

# Make average reference for treatment
ref_obj$celltype.class <- paste0(ref_obj$Celltype_Annotation, ".",
                                 ref_obj$HTO_Classification)
Idents(ref_obj) <- "celltype.class"
ref_mat_new <- AverageExpression(ref_obj, "gene_RNA")
ref_mat_new <- ref_mat_new$gene_RNA

seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat_new,
                                 save_name = "celltype_class",
                                 ADT = ADT,
                                 nfeatures = 1500,
                                 clusters = clusters)

seurat_data <- seurat_res_list$object

seurat_data$RNA_class <- sub(".*\\.", "", seurat_data$RNA_celltype_class)

seurat_data$ADT_class <- sub(".*\\.", "", seurat_data$ADT_celltype_class)

seurat_data$combined_class <- sub(".*\\.", "",
                                  seurat_data$combined_celltype_class)

# Check the classes on the RNA data because this best separates the activation
# levels
pdf(file.path(save_dir, "images", "class_mapping.pdf"))
plots <- plotDimRed(seurat_data, col_by = c("RNA_class",
                                            "ADT_class",
                                            "combined_class"),
                    plot_type = "rna.umap")

print(plots[[1]])
print(plots[[2]])
print(plots[[3]])

dev.off()

## Bulk reference --------------------------------------------------------------
ref_dir <- file.path(base_dir, "files", "celltype_refs", "bulk_RNAseq")

# Ref info
ref_mat <- read.table(file.path(ref_dir,
                             "gene_id_GSE118165_RNA_gene_abundance.txt"))

seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat,
                                 save_name = "celltype_bulk",
                                 ADT = ADT,
                                 nfeatures = 1000,
                                 clusters = clusters)


seurat_data <- seurat_res_list$object

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))