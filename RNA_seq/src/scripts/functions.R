#############
# Functions # 
#############

# Find alignment quality
get_alignment_qual <- function(rate) {
  if (rate >= 85) {
    qual <- "high"
  } else if (rate >= 65) {
    qual <- "okay"
  } else {
    qual <- "low"
  }
  
  qual
}

# Create a PCA Plot
plot_pca <- function(vsd, group_by, color_palette = NULL){
  if(is.null(color_palette)){
    color_palette <- brewer.pal(length(levels(vsd[[group_by]])), "Set1")
    names(color_palette) <- levels(vsd[[group_by]])
  }
  pcaData <- plotPCA(vsd, intgroup = group_by, returnData = T)
  colnames(pcaData)[3] <- "group_by"
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca_plot <- ggplot(data = pcaData,
                     mapping = aes(x = PC1,
                                   y = PC2,
                                   color = group_by)) +
    theme_classic() +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(color = group_by) +
    #coord_fixed() +
    scale_color_manual(values = color_palette)
  return(pca_plot)
}


# Function to determine sample distances
get_distances <- function(vsd){
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  dist_map <- pheatmap(sampleDistMatrix,
                       clustering_distance_rows=sampleDists,
                       clustering_distance_cols=sampleDists,
                       col=colors)
  return(dist_map)
}

# Function to return sig genes for a comparison
get_de <- function(object, column = NULL, var1 = NULL, var2 = NULL,
                   write_csv = F, de_test = NULL, interaction = FALSE,
                   p_value = 0.05, lfc = 0.5, lfc_shrink = TRUE,
                   output_dir = NULL, save_name = NULL){
  if(interaction){
    contrast = de_test
  } else {
    contrast = c(column, var1, var2)
  }
  res <- DESeq2::results(object, contrast = contrast)
  if(lfc_shrink){
    res <- DESeq2::lfcShrink(object, contrast = contrast,
                             res = res, type = "ashr")
  }
  resOrdered <- res[order(res$pvalue), ]
  sigGenes <- resOrdered[!is.na(resOrdered$padj), ]
  sigGenes <- sigGenes[sigGenes$padj < p_value, ]
  sigGenes <- sigGenes[abs(sigGenes$log2FoldChange) > lfc, ]
  # Add columns for gene id and ens id
  sigGenes_colnames <- colnames(sigGenes)
  sigGenes$gene_name <- sub("_ENSMUSG.*", "", rownames(sigGenes))
  sigGenes$ens_id <- sub(".*_ENSMUSG", "ENSMUSG", rownames(sigGenes))
  # Place the new columns at the front
  sigGenes <- sigGenes[ , c("gene_name", "ens_id", sigGenes_colnames)]
  if(write_csv){
    file_name <- paste0(save_name, ".csv")
    write.csv(sigGenes,
              file = file.path(output_dir, "DE_files", file_name),
              row.names = FALSE)
  }
  return(list(DE_genes = sigGenes,
              DE_df = res))
}

# Function to make a heatmap
make_heatmap <- function(dds, vsd, de_genes, treatment, control, group, 
                         print_genenames = FALSE, gene_identifier = "Gene_ID",
                         cluster_cols = FALSE, save_heatmap = TRUE,
                         output_dir = "results", plot_groups = "all",
                         color_test = NULL, save_name = NULL){
  # First make a "col data" data frame. This is just telling
  # pheatmap what you want to use to color the columns.
  df <- as.data.frame(colData(dds)[,c(group)])
  sample_info <- colData(dds)
  colnames(df) <- group
  rownames(df) <- rownames(colData(dds))
  
  # grab genes
  if(isS4(de_genes)){
    genes <- rownames(de_genes)
  } else if (is.character(de_genes)){
    genes <- de_genes
  } else {
    stop("de_genes must be a S4 output from DESeq2 (res$DE_genes) or a list of genes")
  }
  if(is.null(color_test)){
    color_test <- brewer.pal(length(levels(dds[[group]])), "Set1")
    names(color_test) <- levels(dds[[group]])
  }
  coloring <- list(color_test)
  names(coloring) <- group
  heatmap_df <- assay(vsd)[genes,]
  if(!("all" %in% plot_groups)){
    sample_info <- colData(dds)
    sample_plot <- sample_info[sample_info[[group]] %in% plot_groups, ]
    heatmap_df <- heatmap_df[ , rownames(sample_plot)]
  } else {
    sample_plot <- sample_info
  }
  # Heatmaps are always mean cenetered. Meaning the value shown is actually the 
  # expression value of the sample minus the mean. This command centers the
  # dataframe so it can be plotted. When it isn't centered it is very hard to
  # see any trends because all genes have very different expression levels.
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  palOut <- colorRampPalette(blueYellow)(256)
  if(!cluster_cols){
    if(!identical(colnames(heatmap_scale), rownames(sample_plot))){
      heatmap_scale <- heatmap_scale[ , rownames(sample_plot)]
    }
    # Order based on the comparison
    col_order <- c(
      grep(control, sample_plot[[group]]),
      grep(treatment, sample_plot[[group]]),
      grep((paste0(treatment, "|", control)), sample_plot[[group]],
           invert = TRUE)
    )
    heatmap_scale <- heatmap_scale[ , col_order]
  }
  if(print_genenames){
    gene_names <- rownames(heatmap_scale)
    if(gene_identifier == "Gene_ID"){
      gene_names <- sub("_ENSMUSG[0-9]*", "", gene_names)
    } else if (gene_identifier == "ENS_ID"){
      gene_names <- sub("*_ENSMUSG", "ENSMUSG", gene_names)
    }
    rownames(heatmap_scale) <- gene_names
  }
  # Here we make the heatmap
  heatmap <- pheatmap(heatmap_scale, cluster_rows = TRUE,
                      cluster_cols = cluster_cols,
                      show_rownames = print_genenames,
                      show_colnames = TRUE, annotation_col = df,
                      annotation_colors = coloring, color = blueYellow,
                      border_color = NA, clustering_method = "complete",
                      silent = TRUE)
  
  if(save_heatmap){
    comparison <- paste0(save_name, ".pdf")
    pdf(file.path(output_dir, comparison), width = 10, height = 10)
    print(heatmap)
    dev.off()
  }
  
  return(heatmap)
}

# Run gprofiler separate positive and negative
run_gprofiler <- function(gene_table, pos_name, neg_name,
                          custom_bg = FALSE, 
                          correction_method = "gSCS",
                          exclude_iea = FALSE,
                          save_dir = NULL,
                          plot_dir = NULL,
                          evcodes = FALSE){
  pos_sig_table <- gene_table[gene_table$log2FoldChange > 0, ]
  neg_sig_table <- gene_table[gene_table$log2FoldChange < 0, ]
  
  # Here we order by the log fold change so that "ordered_query" can be true
  pos_sig_table <- pos_sig_table[order(pos_sig_table$log2FoldChange),]
  neg_sig_table <- neg_sig_table[order(neg_sig_table$log2FoldChange),]
  pos_sig_genes <- rev(pos_sig_table$ens_id)
  neg_sig_genes <- neg_sig_table$ens_id
  pos_gost <- gost(query = pos_sig_genes,
                   organism = "mmusculus",
                   ordered_query = TRUE,
                   exclude_iea = exclude_iea,
                   user_threshold = 0.05,
                   correction_method = correction_method,
                   custom_bg = custom_bg,
                   evcodes = evcodes)
  
  neg_gost <- gost(query = neg_sig_genes,
                   organism = "mmusculus",
                   ordered_query = TRUE,
                   exclude_iea = exclude_iea,
                   user_threshold = 0.05,
                   correction_method = correction_method,
                   custom_bg = custom_bg,
                   evcodes = evcodes)
  if(evcodes){
    sig_pos_res <- pos_gost$result[pos_gost$result$significant == TRUE, ] %>%
      dplyr::select(!c(intersection, evidence_codes))
    sig_neg_res <- neg_gost$result[neg_gost$result$significant == TRUE, ] %>%
      dplyr::select(!c(intersection, evidence_codes))
  } else {
    sig_pos_res <- pos_gost$result[pos_gost$result$significant == TRUE, ]
    sig_neg_res <- neg_gost$result[neg_gost$result$significant == TRUE, ]
  }
  
  
  if (!is.null(save_dir)){
    pos_save <- paste0(pos_name, "_from_", pos_name, "_vs_", neg_name, ".csv")
    neg_save <- paste0(neg_name, "_from_", pos_name, "_vs_", neg_name, ".csv")
    pos_save_rds <- paste0(pos_name, "_from_", pos_name, "_vs_",
                           neg_name, ".rds")
    neg_save_rds <- paste0(neg_name, "_from_", pos_name, "_vs_",
                           neg_name, ".rds")
    saveRDS(pos_gost, file = file.path(save_dir,pos_save_rds))
    saveRDS(neg_gost, file = file.path(save_dir, neg_save_rds))
    if(nrow(sig_pos_res) > 1){
      sig_pos_res_csv <- apply(sig_pos_res, 2, as.character)
      write.csv(sig_pos_res_csv, file = file.path(save_dir, pos_save))
    }
    if(nrow(sig_neg_res) > 1){
      sig_neg_res_csv <- apply(sig_neg_res, 2, as.character)
      write.csv(sig_neg_res_csv, file = file.path(save_dir, neg_save))
    }
  }
  # These make the plots
  plots <- list()
  plots$C_GOBP <- gost_plots(sig_pos_res, "GO:BP", pos_name)
  plots$C_GOMF <- gost_plots(sig_pos_res, "GO:MF", pos_name)
  plots$C_GOCC <- gost_plots(sig_pos_res, "GO:CC", pos_name)
  plots$C_KEGG <- gost_plots(sig_pos_res, "KEGG", pos_name)
  plots$C_TF <- gost_plots(sig_pos_res, "TF", pos_name)
  plots$T_GOBP <- gost_plots(sig_neg_res, "GO:BP", neg_name)
  plots$T_GOMF <- gost_plots(sig_neg_res, "GO:MF", neg_name)
  plots$T_GOCC <- gost_plots(sig_neg_res, "GO:CC", neg_name)
  plots$T_KEGG <- gost_plots(sig_neg_res, "KEGG", neg_name)
  plots$T_TF <- gost_plots(sig_neg_res, "TF", neg_name)
  
  if (!is.null(plot_dir)){
    save_name <- paste0(pos_name, "_vs_", neg_name, "_separate.pdf")
    
    # This opens up a pdf file. We will save many images into this file
    pdf(file.path(plot_dir, save_name))
    
    plots_list <- lapply(plots, function(x){
      if (!is.null(x)){
        plot(x)
      }
    })
    
    # This closes the pdf
    dev.off()
    
    
  }
  return_list <- list(pos_gost, neg_gost)
  names(return_list) <- c(pos_name, neg_name)
  return_list$plots <- plots
  return(return_list)
  
}

# Run gprofiler on all DE genes
run_gprofiler_all <- function(gene_table, pos_name, neg_name,
                              custom_bg = FALSE, 
                              correction_method = "g_SCS",
                              exclude_iea = FALSE,
                              save_dir = NULL,
                              plot_dir = NULL,
                              save_name = NULL,
                              evcodes = FALSE){
  
  # Because they are now all combined, we don't want to order
  sig_genes <- rev(gene_table$ens_id)
  gost_res <- gost(query = sig_genes,
                   organism = "mmusculus",
                   ordered_query = FALSE,
                   exclude_iea = exclude_iea,
                   user_threshold = 0.05,
                   correction_method = correction_method,
                   custom_bg = custom_bg,
                   evcodes = evcodes)
  
  if(evcodes){
    sig_res <- gost_res$result[gost_res$result$significant == TRUE, ] %>%
      dplyr::select(!c(intersection, evidence_codes))    
  } else {
    sig_res <- gost_res$result[gost_res$result$significant == TRUE, ]
  }
  
  
  if (!is.null(save_dir)){
    save_name_csv <- paste0(save_name, ".csv")
    save_rds <- paste0(save_name, ".rds")
    saveRDS(gost_res, file = file.path(save_dir, save_rds))
    if(nrow(sig_res) > 1){
      sig_res_csv <- apply(sig_res, 2, as.character)
      write.csv(sig_res_csv, file = file.path(save_dir, save_name_csv))
    }
  }
  plot_name <- paste0(pos_name, "_vs_", neg_name)
  return_list <- list(gost_res)
  names(return_list) <- c(plot_name)
  # These make the plots
  plots <- list()
  plots$GOBP <- gost_plots(sig_res, "GO:BP", plot_name)
  plots$GOMF <- gost_plots(sig_res, "GO:MF", plot_name)
  plots$GOCC <- gost_plots(sig_res, "GO:CC", plot_name)
  plots$KEGG <- gost_plots(sig_res, "KEGG", plot_name)
  plots$TF <- gost_plots(sig_res, "TF", plot_name)
  if (!is.null(plot_dir)){
    save_name <- paste0(pos_name, "_vs_", neg_name, ".pdf")
    
    # This opens up a pdf file. We will save many images into this file
    pdf(file.path(plot_dir, save_name))
    
    
    plots_list <- lapply(plots, function(x){
      if (!is.null(x)){
        plot(x)
      }
    })
    
    # This closes the pdf
    dev.off()
    return_list$plots <- plots_list
  }
  return(return_list)
}

gost_plots <- function(results_table, source, title){
  source_results <- results_table[grep(source, results_table$source), ]
  if (nrow(source_results) > 0) {
    source_results <- source_results[order(source_results$precision), ]
    if(nrow(source_results) > 40){
      source_results <- source_results[1:40, ]
    }
    source_results$term_name <- factor(source_results$term_name, levels = 
                                         unique(source_results$term_name))
    source_results$log_padj <- -log10(source_results$p_value)
    go_plot <-ggplot(source_results, aes(x = precision,
                                         y = term_name,
                                         color = log_padj,
                                         size = intersection_size)) +
      geom_point() +
      theme_classic() +
      scale_size(name = "Intersection",
                 range(c(1,max(source_results$intersection_size)))) +
      viridis::scale_color_viridis() +
      theme(text = ggplot2::element_text(size = 10)) +
      ggtitle(paste0(source, ": ", title)) +
      labs(color = "-log10(p-value)") +
      xlab("Precision (proportion of genes)") +
      ylab("Term")
    
    
    
    return(go_plot)
  }
}

# Make gene plots
plot_genes <- function(gene_id, gene_id_list, deseq_obj,
                       intgroup, plot_ggplot = TRUE,
                       color = NULL, return_data = TRUE,
                       print = TRUE, save_path = NULL,
                       normalized = FALSE, vst = NULL){
  if(normalized & is.null(vst)){
    stop(paste0("If wanted normalized values, you must provide a vst object ",
                "where you performed limma::removeBatchEffect"))
  }
  
  # Locate the index of the gene of interest
  index <- grep(paste0("^", gene_id, "$"), gene_id_list)
  if(length(index) == 0){
    print(paste0(gene_id, " not in deseq object"))
    return(NULL)
  } else if(length(index == 1)) {
    counts_plot <- make_plots(index = index, deseq_obj = deseq_obj,
                              intgroup = intgroup, plot_ggplot = plot_ggplot,
                              color = color, return_data = return_data,
                              print = print, save_path = save_path,
                              normalized = normalized, vst = vst)
    return(counts_plot)
  } else {
    # This is in case a gene has two ensembl values
    plot_list <- lapply(index, function(x) make_plots(index = x,
                                                      deseq_obj = deseq_obj,
                                                      intgroup = intgroup,
                                                      plot_ggplot = plot_ggplot,
                                                      color = color, 
                                                      return_data = return_data,
                                                      print = print,
                                                      save_path = save_path,
                                                      normalized = normalized,
                                                      vst = vst))
    return(plot_list)
  }
}
make_plots <- function(index, deseq_obj, intgroup, plot_ggplot,
                       color, return_data, print, save_path,
                       normalized = FALSE, vst = NULL){
  gene_name <- rownames(deseq_obj)[index]
  if(normalized){
    mat <- assay(vst)
    counts_plot <- t(mat[gene_name, , drop = FALSE]) %>% data.frame()
    colnames(counts_plot) <- "count"
    
    counts_plot$Group <- deseq_obj[[intgroup]]
    
  } else {
    counts_plot <- DESeq2::plotCounts(deseq_obj, gene = gene_name,
                                      intgroup = intgroup,
                                      returnData = TRUE) 
    colnames(counts_plot) <- c("count", "Group")
  }

  if(plot_ggplot){
    if(is.null(color)){
      color <- brewer.pal(length(levels(deseq_obj[[intgroup]])), "Set1")
      names(color) <- levels(deseq_obj[[intgroup]])
    }
    
    ggplot_counts_plot <- ggplot2::ggplot(counts_plot,
                                          ggplot2::aes(x = Group,
                                                       y = count)) +
      ggplot2::geom_point(ggplot2::aes(color = Group), size = 3) +
      ggplot2::scale_color_manual(values = color, name = "Group") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         vjust = 1,
                                                         hjust=1)) +
      ggplot2::ggtitle(gene_name)
    if(print){
      print(ggplot_counts_plot)
    }
    if (!is.null(save_path)){
      file_name <- paste0("counts_plot_", gene_name, ".pdf")
      file_path <- file.path(save_path, file_name)
      ggplot2::ggsave(filename = file_path, plot = ggplot_counts_plot)
    }
    if(return_data){
      return(ggplot_counts_plot)
    }
  } else {
    return(counts_plot)
  }
}
