library(ATACseqQC)
library(Rsamtools)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(ChIPpeakAnno)

bamfile <- snakemake@input[[1]]
#dedup_bamfile <- snakemake@input[[1]]
dedup_bamfile <- snakemake@input[[2]]
print(bamfile)
out_dir <- snakemake@params["output_directory"] # one directory for each sample
out_dir_plots <- snakemake@params["output_plot_dir"]
sample <- snakemake@wildcards["sample"]
tss_output <- snakemake@output$tss_score
seqinformation <- readRDS(snakemake@params$seq_info)
genome <- readRDS(snakemake@params$bsgenome)
txs <- readRDS(snakemake@params$transcripts)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))


# out_dir <- "/beevol/home/wellskri/Analysis/Lori_Sussel/Maria_Hansen/220708_bulk_atac/results/atac_qc"
# bamfile <- "/beevol/home/wellskri/Analysis/Lori_Sussel/Maria_Hansen/220708_bulk_atac/results/bowtie2_cutadapt_trim/KO_2_Aligned.sortedByCoord.out.bam"

bamfile.labels_one <- gsub(".bam", "", basename(bamfile))
bamfile.labels_two <- gsub(".bam", "", basename(dedup_bamfile))

# Library complexity ----------------------------------------------------------
print("### Calculating library complexity ###")

reads_dup <- readsDupFreq(bamfile)

pdf(file.path(out_dir_plots, "complexity.pdf"))

estimateLibComplexity(reads_dup)

dev.off()


# Output: 3 columns
# relative.size - the percent of reads
# values - estimator for the number of unique reads represented at least once in a 
# random sample
# reads - relative size multiplied by the total - which comes from the reds_dup file
# I'm assuming that this is the total number of reads.

# More info - https://rdrr.io/bioc/ATACseqQC/src/R/estimateLibComplexity.R
# More info - https://www.rdocumentation.org/packages/preseqR/versions/4.0.0/topics/ds.rSAC.bootstrap
lib_complexity <- estimateLibComplexity(reads_dup, times = 100,
                                        interpolate.sample.sizes = seq(0.1, 1, by = 0.01))


saveRDS(lib_complexity, file = file.path(out_dir, "lib_complexity.rds"))

# Frag size distribution ------------------------------------------------------
print("### Calculating fragment size distribution ###")

pdf(file.path(out_dir_plots, "frag_size.pdf"))

fragSize <- fragSizeDist(bamfile, bamfile.labels_one)

dev.off()

print("### Calculating fragment size distribution ###")

pdf(file.path(out_dir_plots, "dedup_frag_size.pdf"))

fragSize <- fragSizeDist(dedup_bamfile, bamfile.labels_two)

dev.off()

# Nucleosome positioning ------------------------------------------------------

print("### Calculating nucleosome positioning ###")

# Adjust read start sites to account for transposase
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
bamTop100 <- scanBam(BamFile(dedup_bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags


seqlev <- "1" ## subsample data for quick run
#seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Change to match our BAM
# This line frequently fails, give it three tries
attempt_count <- 0
max_attempts <- 3

while (attempt_count < max_attempts) {
  attempt_count <- attempt_count + 1
  
  tryCatch({
    # Code that attempts to download the file
    # Replace the following line with your actual code for downloading the file
    seqlevelsStyle(seqinformation) <- "NCBI"
    
    # If the file download is successful, break out of the loop
    break
  }, error = function(e) {
    # Print the error message (optional)
    cat("Error:", conditionMessage(e), "\n")
    
    # Wait for some time before retrying (optional)
    Sys.sleep(5)
  })
}

# If the maximum number of attempts is reached, exit the program
if (attempt_count == max_attempts) {
  stop("Failed to update the seq levels after multiple attempts.")
}
which <- as(seqinformation[seqlev], "GRanges")
#which <- as(seqinformation, "GRanges")
gal <- readBamFile(dedup_bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(out_dir, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)

#gal1 <- readBamFile(shiftedBamfile, bigFile=TRUE)

# Promotor transcript body score ----------------------------------------------

print("### Calculating promotor transcript body score ###")

# PT score is calculated as the coverage of promoter divided by the coverage of its
# transcript body. PT score will show if the signal is enriched in promoters.

#txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
seqlevelsStyle(txs) <- "NCBI"

pt <- PTscore(gal1, txs)

pt_plot <- data.frame(pt)

pdf(file.path(out_dir_plots, "promotor_transcript_body.pdf"))

plot1 <- ggplot2::ggplot(pt_plot, ggplot2::aes(x = log2meanCoverage,
	                                           y = PT_score)) +
  ggplot2::geom_point() +
  ggplot2::xlab("log2 mean coverage") +
  ggplot2::ylab("Promoter vs Transcript")

print(plot1)

dev.off()

# Nucleosome free region score ------------------------------------------------

print("### Calculating nucleosome free region score ###")

# NFR score is a ratio between cut signal adjacent to TSS and that flanking the
# corresponding TSS. Each TSS window of 400 bp is first divided into 3
# sub-regions: the most upstream 150 bp (n1), the most downstream of 150 bp (n2),
# and the middle 100 bp (nf). Then the number of fragments with 5’ ends
# overlapping each region are calculated for each TSS. The NFR score for each TSS
# is calculated as NFR-score = log2(nf) - log2((n1+n2)/2). A plot can be generated
# with the NFR scores as Y-axis and the average signals of 400 bp window as X-axis,
# very like a MA plot for gene expression data.

nfr <- NFRscore(gal1, txs)

nfr_plot <- data.frame(nfr)

pdf(file.path(out_dir_plots, "nucleosome_free_regions.pdf"))

plot2 <- ggplot2::ggplot(nfr_plot, ggplot2::aes(x = log2meanCoverage,
	                                            y = NFR_score)) +
  ggplot2::geom_point() +
  ggplot2::xlab("log2 mean coverage") +
  ggplot2::ylab("Nucleosome Free Regions score") +
  ggplot2::ggtitle("NFRScore for 200bp flanking TSSs") +
  ggplot2::xlim(c(-10, 0)) +
  ggplot2::ylim(c(-5, 5))

print(plot2)

dev.off()

# TSS enrichment --------------------------------------------------------------

print("### Calculating TSS enrichment ###")

tsse <- TSSEscore(gal1, txs)
tss_score <- tsse$TSSEscore

print(tss_score)
print(tss_output)

tss_score_df <- data.frame("sample" = sample, "score" = tss_score)

print(tss_score_df)

write.table(tss_score_df, file = tss_output, sep = ",")

tsse_df <- data.frame(score = tsse$values, position = 100 * (-9.5:9.5))

pdf(file.path(out_dir_plots, "transcription_start_site_enrichment.pdf"))

plot3 <- ggplot2::ggplot(tsse_df, ggplot2::aes(x = position, y = score)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::xlab("distance to TSS") +
  ggplot2::ylab("aggregate TSS score")

print(plot3)

dev.off()



# Split reads -----------------------------------------------------------------
print("### Splitting reads ###")

# The shifted reads will be split into different bins, namely nucleosome free,
# mononucleosome, dinucleosome, and trinucleosome. Shifted reads that do not 
# fit into any of the above bins will be discarded. Splitting reads is a 
# time-consuming step because we are using random forest to classify the
# fragments based on fragment length, GC content and conservation scores3.

# By default, we assign the top 10% of short reads (reads below 100_bp) as
# nucleosome-free regions and the top 10% of intermediate length reads as 
# (reads between 180 and 247 bp) mononucleosome. This serves as the training 
# set to classify the rest of the fragments using random forest. The number of
# the tree will be set to 2 times of square root of the length of the training set.

## run program for chromosome 1 only
txs <- txs[seqnames(txs) %in% "1"]
#genome <- Mmusculus
## split the reads into NucleosomeFree, mononucleosome, 
## dinucleosome and trinucleosome.
## and save the binned alignments into bam files.
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = out_dir)

# Heatmap and coverage curve for nucleosome positions -------------------------
print("### Making heatmap for nucleosome positions ###")
# By averaging the signal across all active TSSs, we should observe that
# nucleosome-free fragments are enriched at the TSSs, whereas the 
# nucleosome-bound fragments should be enriched both upstream and downstream of 
# the active TSSs and display characteristic phasing of upstream and downstream 
# nucleosomes. Because ATAC-seq reads are concentrated at regions of open 
# chromatin, users should see a strong nucleosome signal at the +1 nucleosome, 
# but the signal decreases at the +2, +3 and +4 nucleosomes.

bamfiles <- file.path(out_dir,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))

pdf(file.path(out_dir_plots, "cumulative_percentage.pdf"))

## Plot the cumulative percentage of tag allocation in nucleosome-free 
## and mononucleosome bam files.
cumulativePercentage(bamfiles[1:2], as(seqinformation["1"], "GRanges"))

dev.off()

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

## estimate the library size for normalization
librarySize <- estLibSize(bamfiles)

print("library size")
print(librarySize)

## calculate the signals around TSSs.
print("### Calculating signal around TSS ###")
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)

## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

pdf(file.path(out_dir_plots, "nucleosome_heatmap.pdf"))

#plot heatmap
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)

dev.off()

## get signals normalized for nucleosome-free and nucleosome-bound regions.
print("### Signal based on nucelosome enrichment ###")
out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")

pdf(file.path(out_dir_plots, "tss_enrichment_by_nucleosome.pdf"))

## rescale the nucleosome-free and nucleosome signals to 0~1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")

dev.off()