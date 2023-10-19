library(DiffBind)
library(tidyverse)
library(GenomicRanges)

hmmratac_addition <- 50


samples <- snakemake@params[["samples"]]
tissue <- snakemake@params[["tissue"]]
bam_files <- snakemake@input[["bam_files"]]
peak_files <- snakemake@input[["peak_files"]]
blacklist <- snakemake@params[["blacklist"]]
output_dir <- snakemake@params[["out_dir"]]
out_obj <- snakemake@output[["db_obj"]]
count_obj <- snakemake@output[["count_obj"]]
wildcards <- snakemake@wildcards

# Set peak caller type
if(wildcards[["peak_caller"]] == "macs2"){
  peakcaller <- "narrow"
} else if (wildcards[["peak_caller"]] == "hmmratac"){
  peakcaller <- "bed"

  # Also need to adjust the narrow summits --------
  new_peak_files <- lapply(peak_files, function(x){
    # read in peak file
    peak_file <- rtracklayer::import(x, format = "bed")

    # Extend the range by 50 nt because the summits file is 1nt range
    peak_file <- peak_file + hmmratac_addition

    save_file <- gsub("_summits", "_summits_extended", x)
    # write file
    rtracklayer::export.bed(peak_file, con = save_file)
    return(save_file)
    })
  peak_files <- unlist(new_peak_files)
} else {
  stop(paste("Unknown peak caller!\nPlease check the specifications of the peak caller and",
             "update the `DiffBind.R` script!"))
}

names(bam_files) <- samples
names(peak_files) <- samples

sample_df <- lapply(samples, function(x){
  sample_id <- x
  condition <- gsub("_[0-9]", "", x)
  replicate <- gsub(".*_", "", x)
  bam_reads <- bam_files[[x]]
  if(!grep(x, bam_reads)){
    exit("bam file isn't from the correct sample")
  }
  peaks <- peak_files[[x]]
  if(!grep(x, peaks)){
    exit("peak file isn't from the correct sample")
  }
  
  tissue_type = tissue[[x]]

  return_df <- data.frame("SampleID" = sample_id,
                          "Tissue" = tissue_type,
                          "Condition" = condition,
                          "Replicate" = replicate,
                          "bamReads" = bam_reads,
                          "Peaks" = peaks,
                          "PeakCaller" = peakcaller)
  
  return(return_df)
})

sample_df <- do.call(rbind, sample_df)

# Create object
diff_bind_obj <- dba(sampleSheet = sample_df)

# downloaded from
# http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/
blacklist <- read.table(blacklist)


blacklist <- GRanges(seqnames = blacklist$V1,
                     ranges = IRanges(start = blacklist$V2,
                                      end = blacklist$V3))

seqlevelsStyle(blacklist) <- "NCBI"


# Remove blacklist -- probably unnecessary for the hmmratac?
diff_bind_obj <- dba.blacklist(diff_bind_obj, greylist = FALSE,
                               blacklist = blacklist)


# Find consensus peaks
olap_rate <- dba.overlap(diff_bind_obj, mode=DBA_OLAP_RATE)

pdf(file.path(output_dir, "consensus_peaks.pdf"))

plot(olap_rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="b")

dev.off()

# Default = peaks in at least 2 samples, change minOverlap to alter
consensus_peaks <- dba.peakset(diff_bind_obj, bRetrieve=TRUE)


# Find counts / peak (even if not a peak in that sample)
atac_counts <- dba.count(diff_bind_obj, peaks = consensus_peaks, score = DBA_SCORE_READS)

# Save count data --> includes FRiP
write.table(x = dba.show(atac_counts), sep = ",",
            file = file.path(output_dir, "sample_info.csv"), quote = FALSE)


# Make correlation plot
pdf(file.path(output_dir, "sample_correlations.pdf"))

plot(atac_counts)

dev.off()

# Make PCA plot, color by condition, label by replicate
pdf(file.path(output_dir, "sample_pca.pdf"))

dba.plotPCA(atac_counts, DBA_CONDITION, label=DBA_REPLICATE)

dev.off()
saveRDS(atac_counts, count_obj)
saveRDS(diff_bind_obj, out_obj)
