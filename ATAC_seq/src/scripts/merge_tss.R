library(tidyverse)

input_files <- snakemake@input
output_file <- snakemake@output[[1]]

all_files <- lapply(intput_files, function(x){
	read.table(x, header = TRUE)
})

all_files = do.call(rbind, all_files)

write.table(all_files, file = output_file, sep = ",")