library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(pathview)
library(openxlsx)
library(gprofiler2)
library(scAnalysisR)
library(here)


pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "sample"
cell_types <- "combined_celltype"
clusters <- "combined_cluster"

gen_id <- "hsa"

path_id_list <- c(NFKB_sig_path = "04064", T_cell_receptor = "04660",
                  cytokine_receptor_interaction = "04060",
                  cell_adhesion = "04514", Hematopoitic_lineage = "04640",
                  Type1_diabetes = "04940", NK_cytoxicity = "04650",
                  th1_th2_differentiation = "04658", TNF_signaling = "04668",
                  immunodeficiency = "05340", Th17_differentiation = "04659")

normalization_method <- "log" # can be SCT or log


if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

source(file.path(base_dir, "src", "scripts", "functions.R"))

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

out_dir <- file.path(save_dir, "images", "pathways")

out_dir_go <- file.path(save_dir, "images", "GSEA")
out_dir_text <- file.path(save_dir, "files", "GSEA")

# Create output directory
out_dir %>%
  dir.create(showWarnings = F)

out_dir_go %>%
  dir.create(showWarnings = F)

out_dir_text %>%
  dir.create(showWarnings = F)


# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))


# Cell type DE -----------------------------------------------------------------

# Find markers
marker_genes <- find_markers(seurat_object = seurat_data,
                             seurat_assay = seurat_assay,
                             test_idents = cell_types)


# Run gost
gost_output <- run_gost(seurat_de = marker_genes,
                        sources = c("GO:BP", "KEGG", "REAC", "TF"))

# Save gost output
save_gost(gost_output,
          save_dir_plots = file.path(out_dir_go, "celltype"),
          save_dir_text = file.path(out_dir_text, "celltype"))



de_to_pathview(seurat_de = marker_genes,
               path_id_list = path_id_list,
               out_dir = file.path(out_dir, "celltype"),
               seurat_assay = seurat_assay,
               test_idents = cell_types,
               gen_id = gen_id)

# cluster DE -------------------------------------------------------------------

# Find markers
marker_genes <- find_markers(seurat_object = seurat_data,
                             seurat_assay = seurat_assay,
                             test_idents = clusters)


# Run gost
gost_output <- run_gost(seurat_de = marker_genes,
                        sources = c("GO:BP", "KEGG", "REAC", "TF"))

# Save gost output
save_gost(gost_output,
          save_dir_plots = file.path(out_dir_go, "cluster"),
          save_dir_text = file.path(out_dir_text, "cluster"))


de_to_pathview(seurat_de = marker_genes,
               path_id_list = path_id_list,
               out_dir = file.path(out_dir, "cluster"),
               seurat_assay = seurat_assay,
               test_idents = clusters,
               gen_id = gen_id)

