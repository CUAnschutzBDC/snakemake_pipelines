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
