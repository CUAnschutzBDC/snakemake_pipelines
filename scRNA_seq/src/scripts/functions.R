create_seurat_object <- function(sample, count_path, ADT = TRUE, hashtag = TRUE,
                                 min_features = 200, min_cells = 3){
  sample_path <- paste0(count_path, "/", sample, "/outs/count/filtered_feature_bc_matrix")
  sample_data <- Read10X(data.dir = sample_path)
  if (ADT){
    sample_object <- CreateSeuratObject(counts = sample_data[["Gene Expression"]],
                                        project = sample, min.cells = min_cells, min.features = min_features)
    if (hashtag){
      protein_data <- sample_data[["Antibody Capture"]]
      hashtag_data <- protein_data[grepl("Hashtag", rownames(protein_data)), ]
      ADT_data <- protein_data[!grepl("Hashtag", rownames(protein_data)), ]
      sample_object[["HTO"]] <- CreateAssayObject(
        counts = hashtag_data[ ,Cells(sample_object)])
      sample_object <- NormalizeData(sample_object, assay = "HTO",
                                     normalization.method = "CLR")
    } else {
      ADT_data <- sample_data[["Antibody_Capture"]]
    }
    sample_object[["ADT"]] <- CreateAssayObject(
      counts = ADT_data[ , Cells(sample_object)])
    sample_object <- NormalizeData(sample_object, assay = "ADT",
                                   normalization.method = "CLR")
  } else {
    sample_object <- CreateSeuratObject(counts = sample_data,
                                        project = sample, min.cells = min_cells,
                                        min.features = min_features)
  }
  return(sample_object)
}

PCA_dimRed <- function(sample_object, assay = "RNA"){
  if(assay == "RNA"){
    DefaultAssay(sample_object) = "RNA"
    sample_object <- RunPCA(sample_object,
                            features = VariableFeatures(object = sample_object))
  } else if(assay == "ADT"){
    DefaultAssay(sample_object) <- "ADT"
    # Use all ADTS for dimensional reduction
    VariableFeatures(sample_object) <- rownames(sample_object[["ADT"]])
    sample_object <- ScaleData(sample_object) %>% 
      RunPCA(reduction.name = "apca")
  } else if(assay == "SCT"){
    DefaultAssay(sample_object) = "SCT"
    sample_object <- RunPCA(sample_object,
                            features = VariableFeatures(object = sample_object),
                            reduction.name = "sctpca")
  }
  return(sample_object)
}

plot_PCA <- function(sample_object, HTO = FALSE, assay = "RNA"){
  if(assay == "RNA"){
    reduction <- "pca"
  } else if(assay == "ADT"){
    reduction <- "apca"
  } else if(assay == "SCT"){
    reduction <- "sctpca"
  } else {
    stop("'reduction' must be 'ADT', 'SCT', or 'RNA'")
  }
  plots <- list()
  plots$pca_loadings <- VizDimLoadings(sample_object, dims = 1:2,
                                       reduction = reduction)
  plots$pca_plot <- plotDimRed(sample_object, plot_type = reduction,
                               col_by = "orig.ident")
  plots$mito_plot <- plotDimRed(sample_object, plot_type = reduction,
                                col_by = "percent.mt")
  plots$nfeature_plot <- plotDimRed(sample_object, plot_type = reduction,
                                    col_by = "nFeature_RNA")
  plots$ncount_plot <- plotDimRed(sample_object, plot_type = reduction,
                                  col_by = "nCount_RNA")
  if(HTO){
    plots$hto_pca_plot <- plotDimRed(sample_object, plot_type = reduction,
                                     col_by = "HTO_classification")
  }
  
  if(assay == "RNA"){
    sample_object <- JackStraw(sample_object, num.replicate = 100,
                               reduction = reduction)
    sample_object <- ScoreJackStraw(sample_object, dims = 1:20)
    plots$jackstraw <- JackStrawPlot(sample_object, dims = 1:20)
  }
  plots$elbow <- ElbowPlot(sample_object, reduction = reduction, ndims = 40)
  return(plots)
}

group_cells <- function(sample_object, sample_name = NULL, save_dir = NULL,
                        nPCs = 10, resolution = 0.8, assay = "RNA",
                        HTO = FALSE, ...){
  if(assay == "RNA"){
    DefaultAssay(sample_object) = "RNA"
    save_plot <- paste0(save_dir, "images/rnaUMAP_", sample_name, ".pdf")
    sample_object <- FindNeighbors(sample_object, dims = 1:nPCs,
                                   reduction = "pca")
    sample_object <- FindClusters(sample_object, resolution = resolution)
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs, reduction = "pca",
                             assay = assay, reduction.key = "rnaUMAP_", reduction.name = "rna.umap")
    sample_object[["RNA_cluster"]] <- Idents(sample_object)
    col_by_list <- c("RNA_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = "rna.umap", ...)
  } else if (assay == "ADT"){
    DefaultAssay(sample_object) <- "ADT"
    save_plot <- paste0(save_dir, "images/adtUMAP_", sample_name, ".pdf")
    sample_object <- FindNeighbors(sample_object,
                                   features = rownames(sample_object),
                                   dims = 1:nPCs, reduction = "apca")
    sample_object <- FindClusters(sample_object, resolution = resolution,
                                  graph.name = "ADT_snn")
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs, reduction = "apca",
                             assay = assay, reduction.key = "adtUMAP_", reduction.name = "adt.umap")
    sample_object[["ADT_cluster"]] <- Idents(sample_object)
    col_by_list <- c("ADT_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = "adt.umap", ...)
  } else if(assay == "SCT"){
    DefaultAssay(sample_object) = "SCT"
    save_plot <- paste0(save_dir, "images/rnasctUMAP_", sample_name, ".pdf")
    sample_object <- FindNeighbors(sample_object, reduction = "sctpca",
                                   dims = 1:nPCs)
    sample_object <- FindClusters(sample_object, resolution = resolution)
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs, reduction = "sctpca",
                             assay = assay, reduction.key = "rnasctUMAP_", reduction.name = "rna_sct.umap")
    sample_object[["SCT_cluster"]] <- Idents(sample_object)
    col_by_list <- c("SCT_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = "rna_sct.umap", ...)
  }
  
  return(list(object = sample_object,
              plots = plot_list))
}

plotDimRed <- function(sample_object, col_by, save_plot = NULL,
                       plot_type = "umap",
                       dims_use = NULL, highlight_group = FALSE,
                       group = NULL, meta_data_col = "orig.ident",
                       return_plot = TRUE, ...) {
  plot_list <- lapply(col_by, function(x) {
    plotDimRedSingle(seurat_object = sample_object, col_by = x, plot_type = plot_type,
                     dims_use = dims_use, highlight_group = highlight_group,
                     group = group, meta_data_col = meta_data_col, ...)
  })
  if (!is.null(save_plot)){
    pdf(save_plot)
    print(plot_list)
    dev.off()
  }
  return(plot_list)
}

plotDimRedSingle <- function(seurat_object, col_by, plot_type = "umap",
                             dims_use = NULL, highlight_group = FALSE,
                             group = NULL, meta_data_col = "orig.ident", ...) {
  # Determine where in Seurat object to find variable to color by
  if (col_by == "cluster" | col_by == "Cluster"){
    col_by_data <- as.data.frame(Idents(object = seurat_object))
  }else if (col_by %in% rownames(seurat_object) |
            col_by %in% colnames(seurat_object[[]])){
    col_by_data <- FetchData(object = seurat_object, vars = col_by)
  }else if (col_by %in% rownames(seurat_object[["ADT"]])){
    col_by_data <- FetchData(object = seurat_object, vars = paste0("adt_", col_by))
  }else {
    stop("col_by must be a gene, metric from meta data or 'cluster'")
  }
  
  # Make the name in the data frame the same regardless of what it was originally
  names(col_by_data) <- "colour_metric"
  
  if (is.null(dims_use)){
    dims_use <- c(1,2)
  }
  # Make a data frame based on the cell embeddings from the plot type of choice
  if (plot_type %in% names(seurat_object)){
    plot_coord <- Embeddings(object = seurat_object, reduction = plot_type)
    plot_names <- colnames(plot_coord)
    ndims <- length(plot_names)
    plot_cols <- lapply(dims_use, function(x){
      if (x > ndims) {
        stop("dims_use must be equal to or less than number of dimensions")
      } else {
        plot_col <- plot_names[x]
        return(plot_col)
      }
    })
    plot_cols <- unlist(plot_cols)
    plot_coord <- plot_coord[ , colnames(plot_coord) %in% plot_cols]
    axis_names <- colnames(plot_coord)
    colnames(plot_coord) <- c("dim1", "dim2")
    plot_df <- merge(plot_coord, col_by_data, by = "row.names")
  } else {
    stop("plot type must be a dimensional reduction in Seurat object")
  }
  # Add in group information if highlighting one group.
  if (highlight_group){
    if (is.null(group)){
      stop("if highlight_group is true, group must be a value from the meta_data
            column specified")
    }
    if (!identical(rownames(seurat_object[[]]), rownames(plot_df))) {
      print("must reorder cells")
      plot_df <- plot_df[match(rownames(seurat_object[[]]),
                               rownames(plot_df)), , drop = FALSE]
    }
    plot_df[[meta_data_col]] <- seurat_object[[meta_data_col]]
    if (is.factor(plot_df$all)){
      plot_df$all <- factor(plot_df$all,
                            levels = c("all_samples", levels(plot_df$all)))
    }
    plot_df$all[!(plot_df[[meta_data_col]] %in% group)] <- "all_samples"
    
    # Plot as descrete
    if (!is.numeric(plot_df$colour_metric)){
      return_plot <- groupDiscretePlots(group, plot_df, axis_names = axis_names,
                                        col_by = col_by, ...)
      # Plot as continuouts
    }else{
      return_plot <- groupContinuousPlots(group, plot_df, axis_names = axis_names,
                                          col_by = col_by, ...)
    }
  }
  # Plot as discrete
  if (!is.numeric(plot_df$colour_metric)){
    return_plot <- discretePlots(plot_df, axis_names = axis_names,
                                 col_by = col_by, ...)
    
    # Plot as continuous
  }else{
    return_plot <- continuousPlots(plot_df, axis_names = axis_names,
                                   col_by = col_by, ...)
  }
  return(return_plot)
}



discretePlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                          color = NULL, save_plot = NULL, show_legend = TRUE,
                          size = 0.25){
  base_plot <- ggplot2::ggplot(data = plot_df,
                               ggplot2::aes_(~dim1, ~dim2)) +
    xlab(axis_names[1]) +
    ylab(axis_names[2])
  
  # Add colors based on metric chosen
  base_plot <- base_plot +
    ggplot2::geom_point(ggplot2::aes_(colour = ~colour_metric),
                        show.legend = show_legend, size = size) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  
  nColors <- length(levels(factor(plot_df$colour_metric)))
  
  # Color based on RColorBrewer if own palette isn't chosen
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_manual(values = color, name = col_by)
  }
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}



continuousPlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                            color = NULL, save_plot = NULL, show_legend = TRUE,
                            size = 0.25){
  base_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes_(~dim1, ~dim2)) +
    xlab(axis_names[1]) +
    ylab(axis_names[2])
  
  base_plot <- base_plot +
    ggplot2::geom_point(ggplot2::aes_(colour = ~colour_metric),
                        show.legend = show_legend, size = size)
  
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_viridis_c(option = "magma") +
      labs(color = col_by)
  } else {
    low <- color[1]
    high <- color[2]
    base_plot <- base_plot + 
      ggplot2::scale_color_gradient(low = low, high = high, name = col_by)
  }
    
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}


groupDiscretePlots <- function(group, plot_df, col_by, axis_names = c("dim1", "dim2"),
                               color = NULL, save_plot = NULL, show_legend = TRUE,
                               size = 0.25) {
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1,
                                                           ~dim2)) +
    xlab(axis_names[1]) +
    ylab(axis_names[2])
  
  base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                               ggplot2::aes_(~dim1, ~dim2), 
                                               color = "#DCDCDC",
                                               size = size,
                                               show.legend = FALSE)
  base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                               ggplot2::aes_(~dim1, ~dim2,
                                                             color = ~all),
                                               size = size,
                                               show.legend = show_legend)
  
  base_plot <- base_plot + #ggplot2::theme_classic() + 
    ggplot2::ggtitle(paste(group, collapse = "_")) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2]) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  if (is.null(color)) {
    nColors <- length(levels(factor(plot2$all)))
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot +
      ggplot2::scale_color_manual(values = color, name = col_by)
  }
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}

groupContinuousPlots <- function(group, plot_df, col_by, color = NULL,
                                 limits = NULL, axis_names = c("dim1", "dim2"),
                                 save_plot = NULL, show_legend = TRUE,
                                 size = 0.25) {
  plot_name_comb <- paste(group, collapse = "_")
  if (is.null(color)) {
    low <- "#00AFBB"
    high <- "#FC4E07"
  }
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes_(~dim1, ~dim2)) +
    xlab(axis_names[1]) +
    ylab(axis_names[2])
  
  base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                               ggplot2::aes_(~dim1, ~dim2), 
                                               color = "#DCDCDC",
                                               size = size,
                                               show.legend = FALSE)
  base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                               ggplot2::aes_(~dim1, ~dim2,
                                                             color = ~colour_metric),
                                               size = size,
                                               show.legend = show_legend)
  
  base_plot <- base_plot + #ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(plot_name_comb, " ", col_by)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(is.null(limits)){
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                           name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_gradient(low = low, high = high, 
                                                           name = col_by, limits = limits)
  }
  
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}


featDistPlot <- function(seurat_object, geneset, cell_cycle = FALSE,
                         plot_type = "violin",
                         color = NULL, sep_by = "cluster", save_plot = NULL,
                         nrow = NULL, ncol = NULL){
  geneset <- setNames(geneset, geneset)
  if (plot_type == "jitter") {
    # Make jitter plots colored by cell cycle stage
    if(cell_cycle){
      gene_list_cycle <- lapply(geneset, function(x) jitterPlot(
        seurat_object = seurat_object, y_val = x, x_val = sep_by,
        col_by = "cycle_phase", color = c("black", "red", "purple")))
      
      # Arrange all plots into one figure
      plot_list <- gridExtra::arrangeGrob(grobs = gene_list_cycle,
                                          nrow = length(geneset))
    } else {
      # Make a jitter plot based on expression of each gene given in the gene
      # set color by stage
      
      gene_list_stage <- lapply(geneset, function(x) jitterPlot(
        seurat_object = seurat_object, y_val = x, x_val = sep_by,
        color = color))
      
      # Make a plot consisting of all plots made above
      plot_list <- gridExtra::arrangeGrob(grobs = gene_list_stage,
                                          nrow = length(geneset))
    }
    
  }
  if (plot_type == "violin" || plot_type == "both") {
    if (plot_type == "both"){
      plot_jitter <- TRUE
    } else {
      plot_jitter <- FALSE
    }
    gene_list_stage <- lapply(geneset, function(x) violinPlot(
      seurat_object = seurat_object, y_val = x, x_val = sep_by,
      color = color, plot_jitter = plot_jitter))
    
    plot_list <- gridExtra::arrangeGrob(grobs = gene_list_stage,
                                        nrow = length(geneset))
    
  }
  if (!(is.null(save_plot))){
    ggplot2::ggsave(plot_list, plot = base_plot)
  }
  return(plot_list)
}

jitterPlot <- function(seurat_object, y_val, x_val,
                       col_by = NULL, color = NULL) {
  plot_data <- plotDF(seurat_object, y_val, x_val,
                      col_by)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes_(~x_value, 
                                                               ~y_value,
                                                               color = ~col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_point() + ggplot2::geom_jitter(shape = 16) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank())
  
  
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_color_manual(values =
                                    (grDevices::colorRampPalette(
                                      RColorBrewer::brewer.pal(9, 
                                                               "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_color_manual(values = color, 
                                                         name = col_by)  
  }
  
  return(plot_base)
}  


violinPlot <- function(seurat_object, y_val, x_val,
                       col_by = NULL, color = NULL,
                       plot_jitter = FALSE) {
  plot_data <- plotDF(seurat_object, y_val, x_val,
                      col_by)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes_(~x_value, 
                                                               ~y_value,
                                                               fill = ~col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank())
  
  if (plot_jitter) {
    plot_base <- plot_base + ggplot2::geom_jitter(shape = 16)
  }
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_fill_manual(values =
                                   (grDevices::colorRampPalette(
                                     RColorBrewer::brewer.pal(9, 
                                                              "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_fill_manual(values = color, 
                                                        name = col_by)  
  }
  
  return(plot_base)
}

plotDF <- function(seurat_object, y_val, x_val,
                   col_by = NULL) {
  # Add y_value to a data frame used for plotting. This value can be a gene
  # or a value from meta data like nGene
  # Determine where in Seurat object to find variable to color by
  if (y_val %in% rownames(seurat_object) |
      y_val %in% colnames(seurat_object[[]])){
    plot_data <- FetchData(object = seurat_object, vars = y_val)
  }else if (y_val %in% rownames(seurat_object[["ADT"]])){
    plot_data <- FetchData(object = seurat_object, vars = paste0("adt_", y_val))
  }else {
    stop("y_val must be a gene, metric from meta data")
  }
  # Name the column
  names(plot_data) <- "y_value"
  
  # Add a column contining the x_value. This should be something discrete
  # Like timepoint or cluster
  if (x_val %in% colnames(seurat_object[[]])) {
    # Should be able to fix this. Take out as and do x_val = then don't need names()
    x_plot_data <- as.data.frame(seurat_object[[]][, x_val, drop = FALSE])
    #x_plot_data <- data.frame("x_value" = seurat_object@meta.data[, x_val,
    #                                                             drop = FALSE])
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else if (x_val == "cluster") {
    x_plot_data <- as.data.frame(Idents(seurat_object))
    #x_plot_data <- data.frame("x_value" = seurat_object@ident)
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  
  # Name the appropriate column of the plotting data
  names(plot_data)[2] <- "x_value"
  plot_data <- plot_data[match(colnames(seurat_object),
                               rownames(plot_data)), ]
  
  # Determine how to color the plot. Default is the x_value but can be any
  # discrete value.
  if (is.null(col_by)) {
    plot_data$col_by <- plot_data$x_value
  } else if (col_by %in% colnames(seurat_object[[]])) {
    plot_data$col_by <- seurat_object[[]][ , col_by]
  } else if (col_by == "cluster") {
    plot_data$col_by <- Idents(seurat_object)
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  return(plot_data)
}

plot_heatmap <- function(seurat_object, gene_list, meta_col,
                         colors = NULL, meta_df = NULL, color_list = NULL,
                         max_val = 2.5, min_val = -2.5, cluster_rows = FALSE,
                         cluster_cols = FALSE, average_expression = FALSE,
                         plot_meta_col = TRUE){
  if(average_expression){
    # Find average expression of genes in clusters
    Idents(seurat_object) <- meta_col
    heatmap_df <- AverageExpression(seurat_object, seurat = FALSE,
                                    group.by = "ident")
    heatmap_df <- heatmap_df$RNA
    # Test if the colnames look like integers
    character_vals <- 
      suppressWarnings(all(!is.na(as.numeric(as.character(colnames(heatmap_df))))))
    if(is.null(meta_df)){
      sample_info <- seurat_object[[meta_col]]
      # Add levels
      if(is.null(levels(sample_info[[meta_col]]))){
        sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
      }
      meta_df <- data.frame(levels(sample_info[[meta_col]]))
      colnames(meta_df) <- meta_col
      if(character_vals){
        rownames(meta_df) <- paste0("X", meta_df[[meta_col]])
      } else {
        rownames(meta_df) <- meta_df[[meta_col]]
      }
      if(is.null(colors)){
        colors <- brewer.pal(length(levels(sample_info[[meta_col]])), "Set1")
        names(colors) <- levels(sample_info[[meta_col]])
      } 
      # make a list for the column labeing
      color_list <- list(colors)
      names(color_list) <- meta_col
    }
  } else {
    # Pull out data and subset to genes of interest
    heatmap_df <- GetAssayData(seurat_object, slot = "data")
  }
  heatmap_df <- heatmap_df[rownames(heatmap_df) %in% gene_list, ]
  heatmap_df <- data.frame(heatmap_df)
  
  heatmap_df <- heatmap_df[order(match(rownames(heatmap_df), gene_list)), ]
  
  
  if(is.null(meta_df)){
    # Make a df for the column labeling
    sample_info <- seurat_object[[meta_col]]
    # Add levels
    if(is.null(levels(sample_info[[meta_col]]))){
      sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
    }
    if(is.null(colors)){
      colors <- brewer.pal(length(levels(sample_info[[meta_col]])), "Set1")
      names(colors) <- levels(sample_info[[meta_col]])
    } 
    # make a list for the column labeing
    coloring <- list(colors)
    names(coloring) <- meta_col
  } else {
    # The sample info and color list must be provided
    sample_info <- meta_df
    coloring <- color_list
  }
  
  # Set cluster order
  
  cluster_order <- levels(sample_info[[meta_col]])
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  # Colors for heatmap (from the ArchR package)
  blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  
  if(!cluster_cols){
    sample_info <- sample_info[order(match(sample_info[[meta_col]],
                                           cluster_order)), , drop = FALSE]
    rownames(sample_info) <- sub("-", ".", rownames(sample_info))
    if(!identical(colnames(heatmap_scale), rownames(sample_info))){
      heatmap_scale <- heatmap_scale[ , rownames(sample_info)]
    }
  }
  
  if(!plot_meta_col){
    sample_info[[meta_col]] <- NULL
  }
  
  # This makes the values more even
  heatmap_scale <- ifelse(heatmap_scale > max_val, max_val, heatmap_scale)
  heatmap_scale <- ifelse(heatmap_scale < min_val, min_val, heatmap_scale)
  
  heatmap <- pheatmap(heatmap_scale, cluster_rows = cluster_rows,
                      cluster_cols = cluster_cols,
                      show_rownames = TRUE,
                      show_colnames = FALSE, annotation_col = sample_info,
                      annotation_colors = coloring, color = blueYellow,
                      border_color = NA, clustering_method = "complete",
                      silent = TRUE)
  return(heatmap)
}

stacked_barplots <- function(seurat_object, meta_col, color = NULL,
                             percent = TRUE, split_by = NULL){
  if(!is.null(split_by)){
    meta_data <- Seurat::FetchData(seurat_object, vars = c(meta_col, split_by)) %>%
      dplyr::rename(meta_col = 1, split_by = 2) %>%
      dplyr::group_by(split_by) %>%
      base::table() %>%
      base::data.frame()
    
  } else {
    meta_data <- Seurat::FetchData(seurat_object, vars = meta_col) %>%
      base::table() %>%
      base::data.frame() %>%
      dplyr::rename(meta_col = 1) %>%
      dplyr::mutate(split_by = "group1") 
  }
  meta_data <- meta_data %>%
    dplyr::mutate(percents = Freq/sum(Freq))
  # Add colors if not provided
  if(is.null(color)){
    color <- brewer.pal(length(unique(meta_data$meta_col)), "Set1")
  }
  if(percent){
    bar_plot <- ggplot2::ggplot(meta_data, ggplot2::aes(x = split_by,
                                                        y = percents,
                                                        fill = meta_col)) +
      ggplot2::ylab("Percent")
  } else {
    bar_plot <- ggplot2::ggplot(meta_data, ggplot2::aes(x = split_by,
                                                        y = Freq,
                                                        fill = meta_col)) +
      ggplot2::ylab("Count")
  }
  bar_plot <- bar_plot +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::xlab(split_by) +
    ggplot2::labs(fill = meta_col)
  
  return(bar_plot)
}



hypergeometric_test <- function(seurat_object, gene_list, DE_table,
                                DE_p_cutoff = 0.05, DE_lfc_cutoff = 0.5,
                                correction_method = "fdr"){
  # Pull out gene list
  hypergeometric_list <- lapply(names(gene_list), function(list_name){
    # Pull out one gene list
    gene_list_one <- gene_list[[list_name]]
    
    DE_list <- lapply(unique(DE_table$cluster), function(cluster_name){
      # Pull out one DE test and only sig genes
      DE_one <- DE_table %>%
        dplyr::filter(cluster == cluster_name &
                        p_val_adj < DE_p_cutoff &
                        avg_log2FC > DE_lfc_cutoff)
      
      # Find number of overlaps
      x <- length(intersect(gene_list_one$V1, DE_one$gene))
      
      # Length of gene list that overlaps with gene in object (total possible genes 
      # to see in comparison)
      m <- length(intersect(gene_list_one$V1, rownames(seurat_object)))
      
      # All genes from object not in list
      n <- length(setdiff(rownames(seurat_object), gene_list_one$V1))
      
      # Total number of genes
      total <- nrow(seurat_object)
      
      # Length of the DE list
      k <- length(DE_one$gene)
      
      # Calculated expected number of genes
      expected_num <- (m*k)/total
      
      # Calcluate the representation factor
      representation <- x/expected_num
      
      # Calculate the p_val
      p_val <- dhyper(x, m, n, k)
      
      return_df <- data.frame(cluster = cluster_name,
                              gene_list = list_name,
                              overlaps = x,
                              gene_list_len = m,
                              DE_length = k,
                              background_len = n,
                              total_genes = total,
                              expected_overlap = expected_num,
                              overrepresentation = representation,
                              p_val = p_val)
      
      return(return_df)
    })
    full_return <- do.call(rbind, DE_list)
    return(full_return)
  })
  
  full_hypergeometric <- do.call(rbind, hypergeometric_list)
  
  # Calculate adjusted p values
  full_hypergeometric$p_adj <- p.adjust(full_hypergeometric$p_val,
                                        method = correction_method)
  
  return(full_hypergeometric)
}

plot_hypergeom <- function(hypergeom_output, colors = NULL, meta_df = NULL,
                           color_list = NULL,
                           cluster_rows = FALSE, cluster_cols = FALSE,
                           color_palette = NULL,
                           breaks = FALSE){
  
  hypergeom_output$log_adj_pval <- -log10(hypergeom_output$p_adj)
  
  hypergeom_output_w <- hypergeom_output %>%
    dplyr::select(c(cluster, gene_list, log_adj_pval)) %>%
    tidyr::pivot_wider(names_from = gene_list, values_from = log_adj_pval) %>%
    base::data.frame()
  
  rownames(hypergeom_output_w) <- hypergeom_output_w$cluster
  hypergeom_output_w$cluster <- NULL
  hypergeom_output_w <- t(hypergeom_output_w)
  
  if(is.null(meta_df)){
    # Make a df for the column labeling
    sample_info <- data.frame(cluster = colnames(hypergeom_output_w))
    rownames(sample_info) <- sample_info$cluster
    # Add levels
    if(is.null(levels(sample_info$cluster))){
      sample_info$cluster <- factor(sample_info$cluster)
    }
    if(is.null(colors)){
      colors <- brewer.pal(length(levels(sample_info$cluster)), "Set1")
      names(colors) <- levels(sample_info$cluster)
    } 
    # make a list for the column labeing
    coloring <- list(colors)
    names(coloring) <- "cluster"
  } else {
    # The sample info and color list must be provided
    sample_info <- meta_df
    coloring <- color_list
  }
  
  # Set cluster order
  
  cluster_order <- levels(sample_info$cluster)
  # Colors for heatmap (from the ArchR package)
  if(is.null(color_palette)){
    color_palette <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                       "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  } else if (color_palette == "blueRed") {
    pal <- colorRampPalette(c("blue", "white", "red"))
    color_palette <- pal(30)
  }
  
  
  if(!cluster_cols){
    sample_info <- sample_info[order(match(sample_info$cluster,
                                           cluster_order)), , drop = FALSE]
    if(!identical(colnames(hypergeom_output_w), rownames(sample_info))){
      hypergeom_output_w <- hypergeom_output_w[ , rownames(sample_info)]
    }
  }
  
  if((breaks)){
    quantile_breaks <- function(xs, n = 30) {
      breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
      breaks[!duplicated(breaks)]
    }
    
    breaks <- quantile_breaks(hypergeom_output_w, n = 30)
  } else {
    breaks <- NULL
  }
  
  heatmap <- pheatmap(hypergeom_output_w, cluster_rows = cluster_rows,
                      cluster_cols = cluster_cols,
                      show_rownames = TRUE, 
                      show_colnames = TRUE, annotation_col = sample_info,
                      annotation_colors = coloring, color = color_palette,
                      border_color = NA, clustering_method = "complete",
                      silent = TRUE, breaks = breaks)
  return(heatmap)
}

make_reference_w_plot <- function(marker_genes, cluster_num,
                                  cut_p_val = 0.05, ref_mat,
                                  group = NULL,
                                  meta_df = NULL,
                                  colors = NULL){
  
  marker_genes <- marker_genes %>%
    dplyr::filter(p_val_adj < cut_p_val & cluster == cluster_num) %>% 
    dplyr::select(gene) %>%
    filter(gene %in% rownames(ref_mat))
  
  ref_mat_genes <- ref_mat %>%
    t %>%
    as.data.frame %>%
    mutate(row_mean = base::rowMeans(.))
  
  ref_mat_genes <- log(ref_mat_genes/ref_mat_genes$row_mean) %>%
    dplyr::select(all_of(marker_genes$gene))
  
  ref_mat_genes$row_mean <- NULL
  
  pivot_cols <- colnames(ref_mat_genes)
  
  if(!is.null(group)){
    meta_df <- meta_df[rownames(ref_mat_genes), ]
    ref_mat_genes$group <- meta_df[[group]]
  } else {
    ref_mat_genes$group <- rownames(ref_mat_genes)
  }
  
  ref_mat_genes_p <- ref_mat_genes %>%
    mutate(cell_type = rownames(.)) %>%
    tidyr::pivot_longer(cols = all_of(pivot_cols),
                        names_to = "gene",
                        values_to = "expression")
  
  ref_mat_genes_p <- ref_mat_genes_p[order(ref_mat_genes_p$group),]
  
  ref_mat_genes_p$cell_type <- factor(ref_mat_genes_p$cell_type,
                                      levels = unique(ref_mat_genes_p$cell_type))
  
  if(is.null(colors)){
    colors <- colorRampPalette(RColorBrewer::brewer.pal(9,
                                                        "Set1"))(length(
                                                          unique(ref_mat_genes_p$group)))
  }
  
  w_plot <- ggplot2::ggplot(ref_mat_genes_p, 
                            ggplot2::aes(x = cell_type, y = expression,
                                         fill = group)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0, color = "red") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::ylab(
      "log(gene expression value / average expression value of all genes)") +
    ggplot2::ggtitle(paste0("cluster ", cluster_num,
                            " markers in immgen ref"))
  return(w_plot)
}



# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x, convert = "human_mouse"){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  if(convert == "human_mouse"){
    to_symbol = "mgi_symbol"
    from_symbol = "hgnc_symbol"
    to_mart = mouse
    from_mart = human
  } else if (convert == "mouse_human"){
    to_symbol = "hgnc_symbol"
    from_symbol = "mgi_symbol"
    to_mart = human
    from_mart = mouse
  }
  
  genesV2 = getLDS(attributes = c(from_symbol), filters = from_symbol,
                   values = x , mart = from_mart, attributesL = c(to_symbol),
                   martL = to_mart, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(genesV2)
}

#' Run differential expression followed by pathview. Can also run pathview on 
#' precomputed DE output
#' 
#' This function allows you to run pathview given a matrix of DE output
#' @param path_id_list A named list of KEGG path ids. The ids are the 5 number
#' ids, the name will be associated with the output file
#' @param seurat_object OPTIONAL A seurat object that has already been
#' normalized. Must provide either a seurat object or a seurat de.
#' @param seurat_de OPTIONAL the returned object from find_markers. If not 
#' provided, find markers will be run. Must provide either a seurat object or
#' a seurat de.
#' @param seurat_assay OPTIONAL which assay to pull the data from when
#' performing DE. Default is "RNA"
#' @param out_dir OPTIONAL an output directory for the pathview plots. Default
#' it "pathview"
#' @param pval OPTIONAL the pvalue cutoff for the adjusted pvalue from the DE
#' test. Default is 0.05.
#' @param ... OPTIONAL arguments supplied to find_markers
#' @keywords pathview
#' @import tidyverse
#' @import pathview
#' @export

de_to_pathview <- function(path_id_list, seurat_object = NULL,
                           seurat_de = NULL,
                           out_dir = "pathview",
                           pval = 0.05, ...){
  if(is.null(seurat_de) && is.null(seurat_object)){
    stop("must provide either a seurat object or a pre made de list")
  } else if(is.null(seurat_de)){
    marker_genes_rna <- find_markers(seurat_object, ...)
  } else {
    marker_genes_rna <- seurat_de
  }
  # Change to wide form
  marker_genes_rna_long <- marker_genes_rna %>%
    dplyr::filter(p_val_adj < pval) %>%
    tidyr::pivot_wider(names_from = cluster,
                       values_from = "avg_log2FC",
                       id_cols = "gene") %>%
    tibble::column_to_rownames(var = "gene") %>%
    as.matrix()
  
  # Run pathview on all paths
  invisible(lapply(names(path_id_list), function(x) 
    run_pathview(gene_matrix = marker_genes_rna_long,
                 path_id = path_id_list[x],
                 path_name = x,
                 out_dir = out_dir)
  ))
}

#' Run differential expression 
#' 
#' This function runs differential expression analysis using Seurat's
#' functions. This is called by several functions but can be used on
#' its own as well.
#' @param seurat_object A seurat object that has already been normalized.
#' @param test_idents OPTIONAL The identities from the metadata used to perform
#' DE. If no identities are provided, the identites already set in the object
#' will be used.
#' @param seurat_assay OPTIONAL which assay to pull the data from when
#' performing DE. Default is "RNA"
#' @param group OPTIONAL to be used when running DE between different samples.
#' The group refers to the column contining the sample info. Default is NULL.
#' @param group_ident_1 OPTIONAL The sample to use as the main sample to test
#' against. Defaults to the first value of the factor
#' @param group_ident_2 OPTIONAL which sample to use as the second group
#' if more than 2 samples are in the seurat object.
#' @keywords differential expression
#' @import tidyverse
#' @export

find_markers <- function(seurat_object,
                        test_idents = NULL, seurat_assay = "RNA",
                        group = NULL, group_ident_1 = NULL,
                        group_ident_2 = NULL){
  if(!is.null(test_idents)){
    Idents(seurat_object) <- test_idents
  }
  if(is.null(group)){
    marker_genes_rna <- FindAllMarkers(seurat_object, assay = seurat_assay)
  } else {
    # Find DE genes between samples per cluster
    if(is.null(group_ident_1)){
      if(is.null(levels(seurat_object[[group]][[1]]))){
        seurat_object[[group]][[1]] <- factor(seurat_object[[group]][[1]])
      }
      group_ident_1 <- levels(seurat_object[[group]][[1]])[1]
    }
    marker_gene_list <- lapply(unique(Idents(seurat_object)), function(x){
      genes <- FindMarkers(seurat_object, ident.1 = group_ident_1,
                           idents.2 = group_ident_2, group.by = group,
                           subset.ident = x)
      genes$cluster <- x
      genes$gene <- rownames(genes)
      return(genes)
    })
    marker_genes_rna <- do.call(rbind, marker_gene_list)
  }
  return(marker_genes_rna)
}

#' Run pathview on single cell data
#' 
#' This function allows you to run pathview given a matrix of DE output
#' @param gene_matrix a matrix of log fold change values. Columns are samples
#' and rownames are genes. Can include NA values. Must be in a matrix format.
#' @param path_id the KEGG 5 number id for a pathway.
#' @param path_name OPTIONAL the name of the KEGG pathway. This is used to
#' name the output files. If this isn't set, output files are named based
#' on the path id
#' @param out_dir OPTIONAL the output directory to copy the files into. This
#' directory will be created if it doesn't already exist. Default is "pathways"
#' @param multi_state OPTIONAL if multiple samples are included should they
#' all be plotted on the same output image. Default is FALSE
#' @param gen_id OPTIONAL the species ID. Default is "mmu" for mouse. Set if 
#' you aren't using mouse
#' @keywords pathview
#' @import tidyverse
#' @import pathview
#' @export

run_pathview <- function(gene_matrix, path_id, path_name = NULL,
                         out_dir = "pathways", multi_state = FALSE,
                         gen_id = "mmu"){
  
  if(is.null(path_name)){
    path_name = path_id
  }
  
  # Create output directory
  out_dir %>%
    dir.create(showWarnings = F)
  
  pathview_out <- pathview(gene.data = gene_matrix,
                           species = gen_id,
                           pathway.id = path_id,
                           gene.idtype = "SYMBOL",
                           kegg.dir = out_dir,
                           multi.state = multi_state,
                           match.data = multi_state,
                           low = list(gene = "#225ea8", cpd = "blue"),
                           mid = list(gene = "white", cpd = "gray"),
                           high = list(gene = "#e31a1c", cpd = "yellow"),
                           na.col = "#bdbdbd")
  
  if(multi_state){
    orig_name <- str_c(gene_id, path_id, ".pathview.png")
    new_name <- str_c(out_dir, "/", gen_id, path_name,
                      ".pathview.png")
    file.rename(orig_name, new_name)
  }
  invisible(lapply(colnames(gene_matrix), function(cluster){
    orig_name <- stringr::str_c(gen_id, path_id, ".pathview.", cluster, ".png")
    new_name <- stringr::str_c(out_dir, "/", gen_id, "_", path_name,
                               ".pathview.", cluster, ".png")
    tryCatch(
      {
        file.rename(orig_name, new_name)
      },
      warning = function(cond){
        output_file <- stringr::str_c(gen_id, path_id, ".pathview.png")
        if(grepl("cannot rename file", cond) && file.exists(output_file)){
          file.remove(output_file)
        }
      }
    )
    
  }))
  
}

#' Run differential expression followed by gost. Can also run gost on 
#' precomputed DE output
#' 
#' This function allows you to run gost given a matrix of DE output or seurat
#' object
#' @param path_id_list A named list of KEGG path ids. The ids are the 5 number
#' ids, the name will be associated with the output file
#' @param seurat_object OPTIONAL A seurat object that has already been
#' normalized. Must provide either a seurat object or a seurat de.
#' @param sources OPTIONAL a list of what sources to plot. Can be any sources
#' returned by gost. Default is GO:BP, KEGG, REAC and TF
#' @param plot_colors OPTIONAL colors for the samples. Default is blue to red.
#' #' @param intersection_cutoff OPTIONAL how many intersecting genes there must
#' be to plot. Default is 5.
#' @keywords gost
#' @import tidyverse
#' @export

run_gost <- function(seurat_de = NULL, seurat_object = NULL,
                     sources = c("GO:BP", "KEGG", "REAC", "TF"),
                     plot_colors = c("blue", "red"),
                     intersection_cutoff = 5){
  if(is.null(seurat_de) && is.null(seurat_object)){
    stop("must provide either a seurat object or a pre made de list")
  } else if(is.null(seurat_de)){
    marker_genes_rna <- find_markers(seurat_object, ...)
  } else {
    marker_genes_rna <- seurat_de
  }
  marker_genes_gost <- marker_genes %>%
    dplyr::filter(p_val_adj < pval)
  
  gene_list <- lapply(unique(marker_genes_gost$cluster), function(x){
    marker_genes_short <- dplyr::filter(marker_genes_gost, cluster == x) %>%
      dplyr::arrange(p_val_adj)
    return(marker_genes_short$gene)
  })
  
  names(gene_list) <- unique(marker_genes_gost$cluster)
  
  gost_output <- gprofiler2::gost(query = gene_list,
                                  organism = "mmusculus",
                                  multi_query = FALSE,
                                  ordered_query = FALSE,
                                  user_threshold = pval,
                                  custom_bg = NULL,
                                  correction_method = "fdr")
  # Separate this out into it's own function I think, a function where you
  # only pass it the gost_output and sources.
  go_plots <- make_go_plots(gost_output = gost_output,
                            sources = sources,
                            plot_colors = plot_colors,
                            intersection_cutoff = intersection_cutoff)
  return(list(gost_output = gost_output, go_plots = go_plots))
}

#' Makes go plots of different terms after running gost. Can take many sources
#' and will make plots for all queries in the gost return
#' 
#' This function makes plots based on different go terms after running gost. 
#' This function is primarily meant to be called by other functions.
#' @param gost_output the output from running gost.
#' @param sources OPTIONAL a list of what sources to plot. Can be any sources
#' returned by gost. Default is GO:BP, KEGG, REAC and TF
#' @param plot_colors OPTIONAL colors for the samples. Default is blue to red.
#' #' @param intersection_cutoff OPTIONAL how many intersecting genes there must
#' be to plot. Default is 5.
#' @keywords gost
#' @import tidyverse
#' @export

make_go_plots <- function(gost_output,
                         sources = c("GO:BP", "KEGG", "REAC", "TF"),
                         plot_colors = c("blue", "red"),
                         intersection_cutoff = 5){
  all_plots <- lapply(unique(gost_output$result$query), function(query){
    source_plots <- lapply(sources,
                           function(source) make_go_plot_single(
                             gost_output = gost_output,
                             gost_query = query,
                             gost_source = source,
                             plot_colors = plot_colors,
                             intersection_cutoff = intersection_cutoff)
    )
    names(source_plots) <- sources
    return(source_plots)
  })
  names(all_plots) <- unique(gost_output$result$query)
  return(all_plots)
}

#' Makes go plots of different terms for one DE test
#' 
#' This function makes plots based on different go terms after running gost. 
#' This function is primarily meant to be called by other functions.
#' @param gost_output the output from running gost.
#' @param gost_query the name of the cell type or sample from DE to plot
#' @param gost_source what source to plot. Can be any sources returned by gost
#' @param plot_colors OPTIONAL colors for the samples. Default is blue to red.
#' @param intersection_cutoff OPTIONAL how many intersecting genes there must
#' be to plot. Default is 5.
#' @keywords gost
#' @import tidyverse
#' @export

make_go_plot_single <- function(gost_output, gost_query, gost_source,
                         plot_colors = c("blue", "red"),
                         intersection_cutoff = 5){
  gost_output_one <- gost_output$result %>%
    dplyr::filter(query == gost_query & source == gost_source &
                    intersection_size >= intersection_cutoff) %>%
    dplyr::distinct(term_name, .keep_all = TRUE) %>%
    dplyr::top_n(-20, wt = p_value) %>%
    dplyr::arrange(precision)
  
  if(nrow(gost_output_one) > 20){
    gost_output_one <- gost_output_one[1:20,]
  }
  
  gost_output_one$term_name <- factor(gost_output_one$term_name,
                                      levels = gost_output_one$term_name)
  
  gost_plot <- ggplot2::ggplot(gost_output_one,
                               ggplot2::aes(x = precision,
                                            y = term_name,
                                            color = -log10(p_value),
                                            size = intersection_size)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradientn(
      colors=grDevices::colorRampPalette(plot_colors)(n = 299)) +
    ggplot2::xlab("GeneRatio") +
    ggplot2::ylab(paste0(gost_source, " term"))
  return(gost_plot)
}

#' Save output from running gost
#' 
#' This function saves both plots and text output from running gost
#' @param gost_output the output from running run_gost function.
#' @param save_dir_plots the directory to save the plots to
#' @param save_dir_text OPTIONAL the directory to save the text to. If left
#' NULL, will default to the same directory as the plots
#' @param save_plots OPTIONAL if plots should be saved. Default is TRUE
#' @param save_text OPTIONAL if text should be saved to a csv file. The whole
#' gost output will be saved. Default is TRUE.
#' @param save_excel OPTIONAL if an excel file should be created. A file will
#' be created for each query and a tab will be inserted for each type of output.
#' @keywords gost
#' @import tidyverse
#' @export
save_gost <- function(gost_output, save_dir_plots, save_dir_text = NULL,
                      save_plots = TRUE, save_text = TRUE, save_excel = TRUE){
  # Save plots
  if(save_plots){
    # Create output directory
    save_dir_plots %>%
      dir.create(showWarnings = F)
    gost_plots <- gost_output$go_plots
    invisible(lapply(names(gost_plots), function(x){
      pdf(paste0(save_dir_plots, x, "GSE.pdf"))
      gost_plots[[x]]
      dev.off()
    }))
  }
  # Save to plot directory if text directory isn't provided
  if(is.null(save_dir_text)){
    save_dir_text <- save_dir_plots
  }
  gost_text <- gost_output$gost_output$result
  gost_text$parents <- NULL
  # Save all results to csv
  if(save_text){
    # Create output directory
    save_dir_text %>%
      dir.create(showWarnings = F)
    write.csv(gost_text,
              paste0(save_dir_text, "all_GSE_results.csv"))
  }
  if(save_excel){
    # Create output directory
    save_dir_text %>%
      dir.create(showWarnings = F)
    invisible(lapply(unique(gost_text$query), function(gost_query){
      
      # Create excel wb
      gene_wb <- createWorkbook()
      
      # Write to excel wb
      full_list <- lapply(unique(gost_text$source), function(gost_source){
        new_df <- gost_text %>%
          dplyr::filter(source == gost_source & query == gost_query)
        gost_source_write <- sub(":", "_", gost_source)
        addWorksheet(gene_wb, gost_source_write)
        writeData(gene_wb, gost_source_write, new_df)
      })
      
      ## Save workbook to working directory
      saveWorkbook(gene_wb,
                   file = paste0(save_dir_text, gost_query,
                                 "GSE_results.xlsx"),
                   overwrite = TRUE)
    }))
  }
}


# Takes a legend from a ggplot object
# Copied from https://gist.github.com/crsh/be88be19233f1df4542aca900501f0fb
gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

# Confusion matrix from Archr
confusionMatrix <- function(i = NULL, j = NULL){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

# From SCENIC
gene_filter_K <- function (exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * 
                           ncol(exprMat), minSamples = ncol(exprMat) * 0.01) {
  outFile_genesKept <- NULL
  dbFilePath <- NULL
  if (class(scenicOptions) == "ScenicOptions") {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  }else {
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }
  if (is.null(dbFilePath)) 
    stop("dbFilePath")
  if (is.data.frame(exprMat)) {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", 
                                   "", methods("AUCell_buildRankings")), collapse = ", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    stop("'exprMat' should be one of the following classes: ", 
         supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if (any(table(rownames(exprMat)) > 1)) 
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  nCountsPerGene <- rowSums(exprMat, na.rm = T)
  nCellsPerGene <- rowSums(exprMat > 0, na.rm = T)
  message("Maximum value in the expression matrix: ", max(exprMat, 
                                                          na.rm = T))
  message("Ratio of detected vs non-detected: ", signif(sum(exprMat > 
                                                              0, na.rm = T)/sum(exprMat == 0, na.rm = T), 2))
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))
  message("\nNumber of genes left after applying the following filters (sequential):")
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > 
                                                      minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", 
          minCountsPerGene)
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > 
                                                      minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", 
          minSamples, " cells")
  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath)
  genesInDatabase <- colnames(getRanking(motifRankings))
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% 
                                                               genesInDatabase)]
  message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
  genesKept <- genesLeft_minCells_inDatabases
  if (!is.null(outFile_genesKept)) {
    saveRDS(genesKept, file = outFile_genesKept)
    if (getSettings(scenicOptions, "verbose")) 
      message("Gene list saved in ", outFile_genesKept)
  }
  return(genesKept)
}

pairwise_markers <- function(seurat_object, meta_col,
                             p_val_cutoff = 0.05, ...){
  Idents(seurat_object) <- meta_col
  combinations <- combn(unique(Idents(seurat_object)), m = 2)
  all_de <- lapply(1:ncol(combinations), function(x){
    ident1 <- combinations[1, x]
    ident2 <- combinations[2, x]
    marker_genes <- FindMarkers(seurat_object, only.pos = TRUE,
                                ident.1 = ident1, ident.2 = ident2, ...)
    marker_genes$cluster_up <- ident1
    marker_genes$cluster_down <- ident2
    marker_genes$gene_name <- rownames(marker_genes)
    return(marker_genes)
  })
  all_de <- do.call(rbind, all_de)
  
  all_de <- dplyr::filter(p_val_adj < 0.05)
  
  return(all_de)
  
}

