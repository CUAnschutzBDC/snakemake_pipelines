library(tidyverse)
library(SingleCellExperiment)

input <- snakemake@input
outfile <- snakemake@output[[1]]
out_rda <- snakemake@output[[2]]

names(input) <- snakemake@params$barcode

print("reading matrix")

mat_list <- lapply(names(input), function(x){
	new_matrix <- read.table(input[[x]], row.names = 1, header = TRUE)

	# Make sure the barcode naming is correct
	if(!(grepl(x, input[[x]]))){
		stop("Barcode doesn't match file name!")
	}

	new_matrix <- new_matrix %>%
	  dplyr::mutate(sample = x) %>%
	  tibble::rownames_to_column("gene")

	return(new_matrix)
})

print("combining matrix")

final_mat <- do.call(rbind, mat_list)

rm(mat_list)

write.csv(final_mat, file = outfile)

# Pull out data frame for each column in the results, these will become
# the assays

print("Making se assays")

`%nin%` = Negate(`%in%`)

keep_columns <- colnames(final_mat)[colnames(final_mat) %nin% c("sample", "gene")]

assay_list <- lapply(keep_columns, function(x){
  return_df <- final_mat %>%
    dplyr::select(c("gene", "sample", all_of(x))) %>%
    dplyr::rename(values = all_of(x)) %>%
    tidyr::pivot_wider(names_from = sample, values_from = values) %>%
    tibble::column_to_rownames("gene")
  return(return_df)
})

rm(final_mat)

print("Making se object")

names(assay_list) <- keep_columns

se <- SingleCellExperiment(assays = assay_list)

saveRDS(se, file = out_rda)