library(DiffBind)
library(tidyverse)
library(GenomicRanges)
library(here)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))


atac_counts <- readRDS(here("results/hmmratac_diffbind_cutadapt_trim/count_obj.rds"))
diff_bind_obj <- readRDS(here("results/hmmratac_diffbind_cutadapt_trim/db_obj.rds"))
blacklist <- "/beevol/home/wellskri/Analysis/ref/blacklist/mm10/mm10.blacklist.bed"

hmmratac_addition <- 50

sample_sheet <- diff_bind_obj$samples %>%
  dplyr::filter(!grepl("_1", SampleID))

# Remove the failed sample
diff_bind_new <- dba(sampleSheet = sample_sheet)

blacklist <- read.table(blacklist)


blacklist <- GRanges(seqnames = blacklist$V1,
                     ranges = IRanges(start = blacklist$V2,
                                      end = blacklist$V3))

seqlevelsStyle(blacklist) <- "NCBI"


# Remove blacklist -- probably unnecessary for the hmmratac?
diff_bind_new <- dba.blacklist(diff_bind_new, greylist = FALSE,
                               blacklist = blacklist)


# Find consensus peaks
olap_rate <- dba.overlap(diff_bind_new, mode=DBA_OLAP_RATE)


plot(olap_rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="b")


# Default = peaks in at least 2 samples, change minOverlap to alter
consensus_peaks <- dba.peakset(diff_bind_new, bRetrieve=TRUE)


# Find counts / peak (even if not a peak in that sample)
atac_counts <- dba.count(diff_bind_new, peaks = consensus_peaks, 
                         score = DBA_SCORE_SUMMIT,
                         summits = TRUE)

plot(atac_counts)

dba.plotPCA(atac_counts, DBA_CONDITION, label=DBA_REPLICATE)

# My own PCA -------------------------------------------------------------------

# Pull out counts 
counts <- dba.peakset(DBA = atac_counts, bRetrieve = TRUE,
                      DataType = "DBA_DATA_FRAME")

counts <- counts %>%
  dplyr::mutate(full_name = paste(CHR, START, END, sep = "_")) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("full_name") %>%
  dplyr::select(!c(CHR, START, END))

pca <- prcomp(t(counts))

pc_vals <- pca$x %>%
  data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::mutate(group = gsub("_.*", "", sample),
                replicate = gsub(".*_", "", sample))

ggplot2::ggplot(data = pc_vals, ggplot2::aes(x = PC1, y = PC2,
                                             color = group)) +
  ggplot2::geom_point(size = 3)

ggplot2::ggplot(data = pc_vals, ggplot2::aes(x = PC1, y = PC2,
                                             color = replicate)) +
  ggplot2::geom_point(size = 3)


# Back to pipeline -------------------------------------------------------------

dba.plotMA(atac_counts, bNormalized=FALSE, sub="Non-Normalized",
           contrast=list(Ctrl=atac_counts$masks$Ctrl,
                         KO=atac_counts$masks$KO))

atac_counts <- dba.normalize(atac_counts)
dba.plotMA(atac_counts, sub="Normalized (Default)",
           contrast=list(Ctrl=atac_counts$masks$Ctrl,
                         KO=atac_counts$masks$KO))

atac_counts <- dba.normalize(atac_counts,
                             library= "RiP")
dba.plotMA(atac_counts, sub="Normalized (Default)",
           contrast=list(Ctrl=atac_counts$masks$Ctrl,
                         KO=atac_counts$masks$KO))

# RiP does look better than full


# Without batch
atac_model <- dba.contrast(atac_counts,
                           reorderMeta=list(Condition="KO"),
                           minMembers = 2)

atac_model <- dba.analyze(atac_model)

dba.show(atac_model, bContrasts = TRUE)

#dba.plotPCA(atac_model, DBA_CONDITION, label=DBA_REPLICATE)

dba.plotMA(atac_model, sub="DE testing")

dba.plotVolcano(atac_model)


# With batch -------------------------------------------------------------------

# First plot
mm <- model.matrix(~Replicate + Condition, atac_counts$samples)
KO <- colMeans(mm[atac_counts$samples$Condition == "KO", ])
Ctrl <- colMeans(mm[atac_counts$samples$Condition == "Ctrl", ])

atac_counts$samples <- atac_counts$samples %>%
  dplyr::mutate(new_rep = ifelse(Replicate == 1, 2, Replicate))

# Batch correct
batch_counts <- limma::removeBatchEffect(counts, atac_counts$samples$Replicate,
                                         design = model.matrix(~atac_counts$samples$Condition))


pca <- prcomp(t(batch_counts))

pc_vals <- pca$x %>%
  data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::mutate(group = gsub("_.*", "", sample),
                replicate = gsub(".*_", "", sample))

ggplot2::ggplot(data = pc_vals, ggplot2::aes(x = PC1, y = PC2,
                                             color = group)) +
  ggplot2::geom_point(size = 3)

ggplot2::ggplot(data = pc_vals, ggplot2::aes(x = PC1, y = PC2,
                                             color = replicate)) +
  ggplot2::geom_point(size = 3)


contrast <- KO - Ctrl

atac_model2 <- dba.contrast(atac_counts,
                            design = "~Replicate + Condition",
                            minMembers = 2,
                            contrast = KO - Ctrl)

atac_model2 <- dba.analyze(atac_model2)

dba.show(atac_model2, bContrasts = TRUE)

dba.plotPCA(atac_model2, DBA_CONDITION, label=DBA_REPLICATE)

dba.plotMA(atac_model2, sub="DE testing")

dba.plotVolcano(atac_model2)

# Their way

atac_model3 <- dba.contrast(atac_counts,
                            design = "~Replicate + Condition",
                            minMembers = 2,
                            contrast = c("Condition", "KO", "Ctrl"))

atac_model3 <- dba.analyze(atac_model3)

dba.show(atac_model3, bContrasts = TRUE)

dba.plotPCA(atac_model3, DBA_CONDITION, label=DBA_REPLICATE)

dba.plotMA(atac_model3, sub="DE testing")

dba.plotVolcano(atac_model3)

