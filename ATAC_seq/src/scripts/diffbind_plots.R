library(DiffBind)
library(tidyverse)
library(GenomicRanges)
library(here)

diffbind_dir <- here("results/diffbind_cutadapt_trim")

diffbind_obj <- readRDS(file.path(diffbind_dir, "db_obj.rds"))

count_obj <- readRDS(file.path(diffbind_dir, "count_obj.rds"))

output_dir <- here("results/R_analysis")

# Make correlation plot
pdf(file.path(output_dir, "images/sample_correlations.pdf"))

plot(count_obj)

dev.off()

# Make PCA plot, color by condition, label by replicate
pdf(file.path(output_dir, "images/sample_pca.pdf"))

dba.plotPCA(count_obj, DBA_CONDITION, label=DBA_REPLICATE)

dev.off()



# MA plot not normalized
dba.plotMA(count_obj, bNormalized=FALSE, sub="Non-Normalized",
           contrast=list(KO=count_obj$masks$KO,
                         Ctrl=count_obj$masks$Ctrl))


# MA plot default normalization

count_obj <- dba.normalize(count_obj)
dba.plotMA(count_obj, sub="Normalized (Default)",
           contrast=list(KO=count_obj$masks$KO,
                         Ctrl=count_obj$masks$Ctrl))


# Default contrast
atac_model <- dba.contrast(count_obj,
                          reorderMeta=list(Condition="KO"),
                          minMembers = 2)
