library(Seurat)
library(dsb)

# Load in raw data
sample_path <- "results/Hashtag_STOC/outs/count/raw_feature_bc_matrix"
sample_data <- Read10X(data.dir = sample_path)

# Load in HTOs and demultiplex
protein_data <- sample_data[["Antibody Capture"]]
hashtag_data <- protein_data[grepl("Hashtag", rownames(protein_data)), ]
ADT_data <- protein_data[!grepl("Hashtag", rownames(protein_data)), ]
sample_object <- CreateSeuratObject(counts = hashtag_data, min.genes = 5,
                                    assay = "HTO")

sample_object[["ADT"]] <- CreateAssayObject(counts = ADT_data)

# Demultiplex
sample_object <- HTODemux(sample_object, assay = "HTO", positive.quantile = 0.99)

# This gave WAY too many positive cells so I'm just using the seurat object as the positive
print(table(sample_object$HTO_classification))

Idents(sample_object) <- "HTO_classification.global"

negative_object <- subset(sample_object, idents = "Negative")

positive_object <- readRDS("results/R_analysis/stoc_doublet.rds")
Idents(positive_object) <- "DF.classifications_0.25_0.005_120"
positive_object <- subset(x = positive_object, idents = "Singlet")
Idents(positive_object) <- "HTO_classification.global"
positive_object <- subset(x = positive_object, idents = "Singlet")
#positive_object[["ADT"]] <- CreateAssayObject(
#  counts = ADT_data[ , Cells(positive_object)])

# Make sure there are no overlapping cells
length(intersect(rownames(negative_object), rownames(positive_object)))

neg_adt_matrix = GetAssayData(negative_object, assay = "ADT", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(positive_object, assay = "ADT", slot = 'counts') %>% as.matrix()
# Use this to normalize
isotype_ctls <- c("ADT-rIgG1", "ADT-rIgG2a", "ADT-mIgG1",
                  "ADT-mIgG2a", "ADT-rIgG2b")
dsb_norm_prot <- DSBNormalizeProtein(
  cell_protein_matrix = positive_adt_matrix,
  empty_drop_matrix = neg_adt_matrix,
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = isotype_ctls)

# Add to object
positive_object[["ADT"]] <- CreateAssayObject(data = dsb_norm_prot)

# Save
saveRDS(positive_object, "results/R_analysis/stoc_dsb_normalized.rds")
# OR find empty droplets using tools explained


# double check remaining cells against my actual cells
# Subset to the cells I use before normalizing