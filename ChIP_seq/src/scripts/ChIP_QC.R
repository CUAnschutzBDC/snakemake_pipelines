library(ChIPQC)

base_dir <- 
  "/Users/wellskr/Documents/Analysis/Lori_sussel/David_Lorberbaum/Sussel_RA_chip/"

library(dplyr)

sample_sheet <- "results/quality_sample_sheet_cutadapt.csv"

annotation <- "mm10"

project <- "RA_chip"

blacklist_file <- paste0(base_dir, "files/mm10.blacklist.bed")


directory <- "results/R_analysis/images"
# Read in data
samples <- read.csv(paste0(base_dir, sample_sheet))

samples$bamReads <- paste0(base_dir, "/results", samples$bamReads)
samples$bamControl <- paste0(base_dir, "/results", samples$bamControl)
samples$Peaks <- paste0(base_dir, "/results", samples$Peaks)
samples$PeakCaller <- "narrow"
samples$PeakFormat <- "narrow"
#samples$ScoreCol <- 5
#samples$LowerBetter <- FALSE

# Changed shifts because it appaers the fragment lengths are long
# Create ChIPQC object 
#chipObj <- ChIPQC(samples, annotation = annotation,
#                  chromosomes = NULL,  shifts = 1:600,
#                  blacklist = blacklist_file)

chipObj <- ChIPQC(samples, annotation = annotation,
                  blacklist = blacklist_file)

## Create ChIPQC report
ChIPQCreport(chipObj,
             reportName = paste0("ChIP_QC_", project),
             reportFolder = directory)

# Changed shifts because it appaers the fragment lengths are long
saveRDS(chipObj, paste0(base_dir, "results/R_analysis/objects/chipqcObj.rda"))
