library(Seurat)
library(harmony)
library(LaCroixColoR)
library(cowplot)
source("src/scripts/functions.R")
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

stoc_object <- readRDS("results/R_analysis/stoc_processed.rds")


colors <- lacroix_palette("Coconut", 6, "discrete")

colors <- colors[c(1,2,3,5)]

names(colors) <- c("DC2", "Fb_cycling", "Mast", "mTEC.III.")

# Run harmony on the ADTs
stoc_object <- RunHarmony(object = stoc_object,
                          group.by.vars = "HTO_classification",
                          reduction = "apca",
                          dims.use = 1:16,
                          reduction.save = "aharmony",
                          assay.use = "ADT")

stoc_object <- FindNeighbors(stoc_object, dims = 1:20,
                               reduction = "aharmony")
stoc_object <- FindClusters(stoc_object, resolution = 0.6)
stoc_object <- RunUMAP(stoc_object,
                       metric = "correlation", dims = 1:20,
                       reduction = "aharmony",
                       assay = assay,
                       reduction.key = "aharmonyUMAP_",
                       reduction.name = "aharmony.umap")

plotDimRed(sample_object =  stoc_object,
           save_plot = NULL,
           col_by = "HTO_classification",
           plot_type = "aharmony.umap")

stoc_object <- FindMultiModalNeighbors(
  stoc_object, reduction.list = list("pca", "aharmony"), 
  dims.list = list(1:27, 1:20),
  modality.weight.name = c("hRNA.weight", "hADT.weight"),
  knn.graph.name = "harmonywknn",
  snn.graph.name = "harmonywsnn",
  weighted.nn.name = "harmonyweighted.nn"
)

stoc_object <- RunUMAP(stoc_object, nn.name = "harmonyweighted.nn",
                     reduction.name = "wnnharmony.umap",
                     reduction.key = "wnnarmonyUMAP_")
stoc_object <- FindClusters(stoc_object,
                          graph.name = "harmonywsnn",
                          algorithm = 3, resolution = 0.6, verbose = FALSE)

stoc_object[["combined_cluster_harmony"]] <- Idents(stoc_object)

plotDimRed(sample_object = stoc_object,
          save_plot = NULL,
          col_by = "HTO_classification",
          return_plot = TRUE,
          plot_type = "wnnharmony.umap")

violin1 <- featDistPlot(stoc_object, "hRNA.weight",
                        sep_by = "combined_cluster_harmony")
violin2 <- featDistPlot(stoc_object, "hADT.weight",
                        sep_by = "combined_cluster_harmony")
violin3 <- featDistPlot(stoc_object, "hRNA.weight",
                        sep_by = "combined_celltype", color = colors)
violin4 <- featDistPlot(stoc_object, "hADT.weight",
                        sep_by = "combined_celltype", color = colors)
plot_grid(violin1, violin2, violin3, violin4,
          nrow = 2, ncol = 2, align = "hv", axis = "tb")

saveRDS(stoc_object, "results/R_analysis/STOC_harmony.rds")
