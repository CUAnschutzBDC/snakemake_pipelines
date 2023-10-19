library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(here)
library(openxlsx)
library(tidyverse)

#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Functions --------------------------------------------------------------------
analyze_data <- function(sample, peak_granges, peak_path,
                         promoter, save_dir, txdb,
                         de_genes_up, de_genes_down){
  
  tss_title <- paste0("Distribution of ", sample, 
                      " atac peaks relative to TSS")
  
  
  pdf(file.path(save_dir, "coverage_plot.pdf"))
  # Makes a coverage plot across the genome
  print(covplot(peak_granges, weightCol = "score"))
  
  dev.off()
  
  # Gets positions of all promoters
  tagMatrix <- getTagMatrix(peak_granges, windows=promoter)
  
  pdf(file.path(save_dir, "tss_heatmap.pdf"))
  # Heatmap of tss enrichment
  tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
  
  dev.off()
  
  pdf(file.path(save_dir, "tss_enrichment.pdf"))
  # Tss enrichment
  print(plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
  
  dev.off()
  #plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
  
  # Annotates peaks, can change the size of the tss region
  peakAnno <- annotatePeak(peak_path, tssRegion=c(-3000, 3000),
                           TxDb=txdb) 
  
  
  pdf(file.path(save_dir, "peak_enrichment.pdf"))
  # Find enrichment of peaks in different genomic regions
  # Make a handful of useful plots
  print(plotAnnoPie(peakAnno))
  print(plotAnnoBar(peakAnno))
  print(upsetplot(peakAnno))
  
  dev.off()
  
  pdf(file.path(save_dir, "enrichment_up_downstream.pdf"))
  # What percent of binding sites up and downstream of the TSS fall within
  # a certain distance
  print(plotDistToTSS(peakAnno,
                title = tss_title))
  
  dev.off()
  
  # Pull out genes near the peak
  anno_df <- as.data.frame(peakAnno)
  
  gene_mapping <- gene_list$mgi_symbol
  names(gene_mapping) <- gene_list$ensembl_gene_id
  
  anno_df$gene_symbol <- gene_mapping[anno_df$geneId]
  
  # Subset the anno df to only genes that are upregulated when grg is
  # knocked out
  anno_df_up <- anno_df %>%
    dplyr::filter(gene_symbol %in% rownames(de_genes_up))
  
  # Subset the anno df to only genes that are downregulated when grg is 
  # knocked out
  anno_df_down <- anno_df %>%
    dplyr::filter(gene_symbol %in% rownames(de_genes_down))
  
  # Subset the anno df to only genes that aren't changed
  anno_df_none <- anno_df %>%
    dplyr::filter(gene_symbol %notin% rownames(beta_de_genes))
  
  # TSS plot separated
  anno_df_up$.id <- paste0(sample, "_up")
  anno_df_down$.id <- paste0(sample, "_down")
  anno_df_none$.id <- paste0(sample, "_none")
  
  # return this as well
  anno_df_all <- do.call(rbind, list(anno_df_up, anno_df_down, anno_df_none))
  
  pdf(file.path(save_dir, "enrichment_by_DE.pdf"))
  print(ChIPseeker:::plotDistToTSS.data.frame(anno_df_all,
                                        distanceColumn = "distanceToTSS",
                                        categoryColumn = ".id",
                                        title = tss_title))

  dev.off()
  
  return(list("peak_annotation" = peakAnno, 
              "peaks_with_DE" = anno_df_all))
  
}

`%notin%` <- Negate(`%in%`) 

# Read in data -----------------------------------------------------------------

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Make txdb object from the GTF used to do the alignment
txdb <- makeTxDbFromGFF(file.path("/Users/wellskr/Documents/Analysis/references",
                                  "gtf/mouse/Mus_musculus.GRCm38.96.gtf"),
                        organism = "Mus musculus")

macs_dir <- here("results/macs2_cutadapt_trim")

de_dir <- file.path("/Users/wellskr/Documents/Analysis/Lori_sussel/",
                    "Alex_Theis/sussle_220203_mouse_islet/results_tomato",
                    "AT_combined_all/R_analysis/files")

samples <- c("Ctrl_1", "Ctrl_2", "KO_1", "KO_2")

# Gene list to get gene ID from ENS ID
gene_list <- readRDS(here("files/ens_to_geneid.rds"))

# get all narrow peak files
files <- lapply(samples, function(x){
  return(file.path(macs_dir, paste0(x, "_peaks.narrowPeak")))
})

names(files) <- samples

# Make all peak files into peak granges objects
peaks <- lapply(files, function(x){
  peak <- readPeakFile(x)
  colnames(GenomicRanges::elementMetadata(peak)) <- c("name", "score", "strand",
                                                      "signalValue",
                                                      "pvalue_neg_log10",
                                                      "qvalue_neg_log10",
                                                      "peak")
  
  return(peak)
})


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# Download beta gene list from the scRNA-seq
beta_de_genes <- openxlsx::read.xlsx(xlsxFile = file.path(de_dir,
                                                          "psuedobulk_de_alone.xlsx"),
                                     sheet = "beta", rowNames = TRUE)

# Find beta genes that are upregulated when grg is knocked out
beta_genes_up <- beta_de_genes %>%
  dplyr::filter(log2FoldChange > 0)

# Beta genes that are down regulated when grg is knocked out
beta_genes_down <- beta_de_genes %>%
  dplyr::filter(log2FoldChange < 0)


chipseeker_dir <- here("results", "chipseeker")
ifelse(!dir.exists(chipseeker_dir), dir.create(chipseeker_dir), FALSE)

all_data <- lapply(samples, function(x){
  save_dir <- file.path(chipseeker_dir, x)
  ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)
  analyze_data(sample = x, peak_granges <- peaks[[x]],
               peak_path = files[[x]], save_dir = save_dir,
               txdb = txdb, de_genes_up = beta_genes_up,
               de_genes_down = beta_genes_down, promoter = promoter)
})

names(all_data) <- samples

all_peak_anno <- lapply(names(all_data), function(x){
  frequencies <- all_data[[x]]$peak_annotation@annoStat
  frequencies$sample <- x
  return(frequencies)
})

all_peak_anno <- do.call(rbind, all_peak_anno)

all_peak_anno <- all_peak_anno %>%
  dplyr::filter(sample != "KO_1")


all_colors <- grDevices::colorRampPalette(colors =
                                            RColorBrewer::brewer.pal(n = 9,
                                                                     name = "Set1"))

feature_palette <- all_colors(length(levels(all_peak_anno$Feature)))
names(feature_palette) <- levels(all_peak_anno$Feature)

pdf(file.path(chipseeker_dir, "all_peak_enrichment.pdf"))

print(ggplot2::ggplot(all_peak_anno, ggplot2::aes(x = sample, y = Frequency,
                                            fill = Feature)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::scale_fill_manual(values = feature_palette))

dev.off()

all_peaks_with_de <- lapply(all_data, function(x){
  return(x$peaks_with_DE)
})


all_peaks_with_de <- do.call(rbind, all_peaks_with_de)

all_peaks_with_de <- all_peaks_with_de %>%
  dplyr::filter(!grepl("KO_1", .id))


pdf(file.path(chipseeker_dir, "enrichment_by_DE.pdf"))
print(ChIPseeker:::plotDistToTSS.data.frame(all_peaks_with_de,
                                            distanceColumn = "distanceToTSS",
                                            categoryColumn = ".id",
                                            title = tss_title))

dev.off()



# Find number of genes in the promoters

promoter_genes <- lapply(all_data, function(x){
  promoter_peaks <- x$peaks_with_DE %>%
    dplyr::filter(annotation == "Promoter (<=1kb)")
  
  promoter_genes <- unique(promoter_peaks$gene_symbol)
  print(length(promoter_genes))
  return(promoter_genes)
  
})

comparisons <- list(c("Ctrl_1", "Ctrl_2"),
                    c("Ctrl_1", "KO_2"),
                    c("Ctrl_2", "KO_2"),
                    c("Ctrl_1", "Ctrl_2", "KO_2"))

intersecting_genes <- lapply(comparisons, function(x){
  intersections <- Reduce(intersect, promoter_genes[x])
  print(length(intersections))
  return(intersections)
})


# To make my own plots eventually ----------------------------------------------

# anno_df_up_right <- anno_df_up %>%
#   dplyr::mutate(distanceToTSS_abs = abs(distanceToTSS)) %>%
#   dplyr::mutate(Feature =
#                   case_when(distanceToTSS == 0 ~ "0kb",
#                             distanceToTSS > 0 & distanceToTSS < 1000 ~ "0-1kb",
#                             distanceToTSS >= 1000 & distanceToTSS < 3000 ~ "1-3kb",
#                             distanceToTSS >= 3000 & distanceToTSS < 5000 ~ "3-5kb",
#                             distanceToTSS >= 5000 & distanceToTSS < 1000 ~ "5-10kb",
#                             distanceToTSS >= 1000 & distanceToTSS < 100000 ~ "10-100kb",
#                             distanceToTSS >= 100000 ~ ">100kb")) %>%
#   dplyr::group_by(Feature) %>%
#   dplyr::count(name = "count") %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(freq = round(count / sum(count), 2) * 100)
# 
# ggplot2::ggplot(anno_df_up_right, ggplot2::aes())