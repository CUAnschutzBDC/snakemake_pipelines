library(tidyverse)

#samples <- snakemake@params["samples"]
# input_files <- snakemake@input
# output_file <- snakemake@output[[1]]
args <- commandArgs(trailingOnly = TRUE)
output_file <- args[[1]]
input_files <- args[2:length(args)]

samples <- stringr::str_split(input_files, pattern = "/")

samples <- lapply(samples, function(x){
  x <- x[length(x) - 1]
})

names(input_files) <- samples

lib_complexity <- lapply(samples, function(x){
  complexity <- readRDS(input_files[[x]])
  complexity$sample <- x
  return(complexity)
})

lib_complexity <- do.call(rbind, lib_complexity)


# Find the sample with the lowest coverage -------------------------------------

one_vals <- lib_complexity %>%
  dplyr::filter(relative.size == 1) %>%
  dplyr::mutate(min_val = min(values)) %>%
  dplyr::mutate(percent_min = round(min_val / values, 2))

min_value <- one_vals$min_val[1]

# Find the downsampling based on the closest value -----------------------------

best_match <- lapply(samples, function(x){
  print(x)
  sample_complexity <- lib_complexity %>%
    dplyr::filter(sample == x) %>%
    dplyr::mutate(difference = abs(values - min_value)) %>%
    dplyr::filter(difference == min(difference))
  
  return(sample_complexity)
})

best_match <- do.call(rbind, best_match)

write.csv(best_match, output_file)