library(Seurat)
library(dsb)
library(tidyverse)
library(here)
library(scAnalysisR)


sample <- "sample"

HTO <- TRUE

# Set directories
base_dir <- here()

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

# Load in raw data
sample_path <- file.path(base_dir, "results", sample,
                      "outs", "count", "raw_feature_bc_matrix")

sample_data <- Read10X(data.dir = sample_path)

if(HTO){

  # Load in HTOs and demultiplex
  protein_data <- sample_data[["Antibody Capture"]]
  hashtag_data <- protein_data[grepl("Hashtag", rownames(protein_data)), ]
  ADT_data <- protein_data[!grepl("Hashtag", rownames(protein_data)), ]
  sample_object <- CreateSeuratObject(counts = hashtag_data,
                                      assay = "HTO")

  sample_object[["ADT"]] <- CreateAssayObject(counts = ADT_data)

  sample_object <- NormalizeData(sample_object, assay = "HTO",
                                 normalization.method = "CLR")

  sample_object <- NormalizeData(sample_object, assay = "ADT",
                                 normalization.method = "CLR")

  sample_object_negative <- subset(sample_object, nCount_HTO == 0)

  sample_object_test <- subset(sample_object, nCount_HTO > 0)



  # Demultiplex
  sample_object_test <- HTODemux(sample_object_test, assay = "HTO",
                                 positive.quantile = 0.99, verbose = TRUE)

  sample_object_negative$HTO_classification.global <- "Negative"

  # This gave WAY too many positive cells so I'm just using the seurat object as the positive
  print(table(sample_object_test$HTO_classification))

  sample_object <- merge(sample_object_negative, sample_object_test, merge.data = TRUE)

  Idents(sample_object) <- "HTO_classification.global"

  # Create an object of just negative cells
  negative_object <- subset(sample_object, idents = "Negative")

  # Read in my exisitng object as the "positive"
  positive_object <- readRDS(file.path(save_dir, "rda_obj", "seurat_doublet.rds"))

  # Make sure there are no overlapping cells
  length(intersect(rownames(negative_object), rownames(positive_object)))

  neg_adt_matrix <- GetAssayData(negative_object, assay = "ADT",
                                 slot = 'counts') %>% as.matrix()
  positive_adt_matrix <- GetAssayData(positive_object, assay = "ADT",
                                      slot = 'counts') %>% as.matrix()

} else {
  # Read in my exisitng object as the "positive"
  starting_object <- readRDS(file.path(save_dir, "rda_obj", "seurat_unprocessed.rds"))

  # define a vector of cell-containing barcodes and remove them from unfiltered data 
  stained_cells <- colnames(starting_object)
  background <- setdiff(colnames(sample_data$`Gene Expression`), stained_cells)

  # split the data into separate matrices per assay 
  prot <- sample_data$`Antibody Capture`
  rna <- sample_data$`Gene Expression`

  # create metadata of droplet QC stats used in standard scRNAseq processing
  rna_size <- log10(Matrix::colSums(rna))
  prot_size <- log10(Matrix::colSums(prot))
  ngene <- Matrix::colSums(rna > 0)
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  propmt <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  md <- as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
  md$bc <- rownames(md)
  md$droplet_class <- ifelse(test = md$bc %in% stained_cells, yes = 'cell',
                            no = 'background')

  # filter barcodes to only include those with data for both assays 
  md <- md %>% dplyr::filter(rna_size > 0 & prot_size > 0 )

  ggplot(md, aes(x = log10(ngene), y = prot_size )) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~droplet_class) 

  # define a vector of background droplet barcodes based on protein library size and mRNA content
  background_drops <- md[md$prot_size > 1.5 &
                          md$prot_size < 3 & md$ngene < 100, ]$bc

  neg_adt_matrix <- as.matrix(prot[ , background_drops])

  # Get positive cells
  positive_object <- readRDS(paste0(save_dir, "rda_obj/seurat_doublet.rds"))

  positive_cells <- colnames(positive_object)
  positive_adt_matrix <- as.matrix(prot[ , positive_cells])

}

# calculate quantiles of the raw protein matrix 
d1 <- data.frame(pmax = apply(positive_adt_matrix, 1, max)) %>% 
  rownames_to_column('prot') %>% arrange(pmax)

# CHECK D1 AND DETERMINE IF ANYTHING NEEDS TO BE REMOVED!!!
# remove non staining CD137_A0355 protein 
#prot_names <- rownames(positive_adt_matrix)
#positive_adt_matrix <- positive_adt_matrix[!prot_names == 'CD137_A0355', ]
#neg_adt_matrix <- neg_adt_matrix[!prot_names == 'CD137_A0355', ]

# Run the normalization
dsb_norm_prot <- DSBNormalizeProtein(
  cell_protein_matrix = positive_adt_matrix,
  empty_drop_matrix = neg_adt_matrix,
  denoise.counts = TRUE,
  use.isotype.control = FALSE)

# Add to object
positive_object[["ADT"]] <- CreateAssayObject(data = dsb_norm_prot)

positive_object[["RAW_ADT"]] <- CreateAssayObject(data = positive_adt_matrix)


# Save
saveRDS(positive_object, file.path(save_dir, "rda_obj", "seurat_doublet.rds"))
