# Single Cell Long Read
Snakemake pipeline to analyze single cell long read data

## To run

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install snakemake -c bioconda -c conda-forge
```

3. Update the config file (config.yaml) 
>* RAW_DATA: Location of the raw data 
>* SAMPLES: the list of samples you want to test. This is the name that will be in the output files.
>* RESULTS: Path to the output directory
>* GENOME: Path to the cellranger genome reference
>* MAX_JOBS: The maximum number of jobs that can be submitted by cell ranger at a time
>* SICELORE_PATH: The path to the [sicelore package](https://github.com/ucagenomix/sicelore)
>* MINION_QC_PATH: The path to [minion qc](https://github.com/roblanf/minion_qc)
>* TENX_DIR: Path to the directory containing output from cellranger
>* CHR_LIST: The list of chromosomes of interest
>* CUTOFF: The fraction of the full transcript length to determine if a transcript should be output as "full length" in the final bam file (step 16).

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`

6. I highly recommend looking at the csv files that are generated and passed to cell ranger to ensure that the correct fastq files have been detected for each sample.



## Steps

1. Use scilore to parse the 10x files and generate barcode, umi, gene sets
2. Run minion QC for nanopore runs for each sample
3. Run read scanner from scilore to throw away reads without polyA tails and 10x adaptors
4. Concatanate all remaining fastq files
5. Split these files for quicker processing
6. Count tags from the scilore pipeline in all reads and in reads thrown out to get an idea for why reads are thrown out. Custom python script that parses the bam files
7. Find the read lengths associated with each category of read to get a better sense of the strucuture of the data. Custom python script that parses the bam files
8. Map reads with minimap 2, tag bams with barcode, gene, and UMI info using scilore
9. Merge all files
10. Create a consensus sequence using UMIs
11. Repeat mapping and tagging on consensus
12. Create count matrix using original, non-deduplicated reads - this is for showing the "counts" per cell
13. Create count matrix using deduplicated reads - for all other downsteram analysis
14. Count using my own strategy, custom python scripts and minimap2
  * Take information from the tagged bams and make fastq files with barcode in the name
  * Align these to the transcriptome allowing for multimappers with identical scores
  * Count only reads that aligned uniquely, these are likely from only one isoform, use read name info to determine the cell that it came from
15. Count using salmon, custom python scripts, minimap2, and salmon
  * Take minimap transcriptome output and split bam files to one bam file per cell (using the read name)
  * Pass each of these bam files to salmon
  * Merge salmon outputs into one matrix per sample
16. Generate files for GVIZ
  * Remove multiple alignments
  * Subset to only "long reads", reads defined as at least x percent of the full length reference isoform
  * These can be used for plotting alignments

## Snakefiles

* `parse_illumina.snake`
  * `parse_illumina`
    * Step 1 - Use scilore to parse the 10x files and generate barcode, umi, gene sets
* `minion_qc.snake`
  * `minion_qc`
    * Step 2 - Run minion QC for nanopore runs for each sample
* `np_read_scanner.snake`
  * `read_scanner`
    * Step 3 - Run read scanner from scilore to throw away reads without polyA tails and 10x adaptors
  * `cat_files`
    * Step 4 - Concatanate all remaining fastq files
  * `split_files`
    * Step 5 - Split these files for quicker processing
  * `minimap.snake`
    * `make_minimap_ref`
      * Make reference for minimap using the cellranger reference
    * `run_minimap`
      * Step 8 - Map reads with minimap 2 against genome reference
    * `minimap2_consensus`
      * Step 11 - Repeat mapping on consensus
    * `make_transcriptome_fa`
      * Step 14 - Make fasta file based on the genome fasta and gtf files
    * `make_minimap_transcript_ref`
      * Step 14 - Make transcriptome reference for minimap
    * `run_mimimap_transcriptome`
      * Step 14 - Map this fastq file against the transcriptome reference using minimap2. Allow it to output minimappers with the same quality. This will allow for equivalent mapping between isoforms when they can't be distinguished.
 * `tag_bams.snake`
   * `tag_genes`
     * Step 8 - Tag minimap output bam with genes using scilore
   * `tag_quality`
     * Step 8 - Tag minimap output bam with quality using scilore
   * `tenx_barcode`
     * Step 8 - Tag minimap output bam with barcode from 10x
   * `count_tags`
     * Step 6 - Count tags from the scilore pipeline in all reads and in reads thrown out to get an idea for why reads are thrown out.
     * output - `{results}/tenx_barcode/{sample}/counted_barcodes.tsv`
   * `find_lengths`
     * Step 7 - Find the read lengths associated with each category of read to get a better sense of the strucuture of the data.
     * output - `{results}/read_lengths/{sample}`
   * `merge_bams`
     * Step 9 - merge all bam files for one sample
   * `bam_all_for_gviz`
     * Step 16 - Subset bam file to only unique reads so gviz will work
   * `tag_genes_consensus`
     * Step 11 - Repeat gene tagging on bam file from consensus minimap
   * `tag_UMI_BC_consensus`
     * Step 11 - Repeat UMI and barcode tagging on bam file from consensus minimap
 * `consensus_sequences.snake`
   * `consensus`
     * Step 10 - Create a consensus sequence using UMIs from scilore
 * `create_matrix.snake`
   * `create_matrix`
     * Step 13 - Create count matrix using deduplicated reads - for all other downsteram analysis
     * output - `{results}/IsoformMatrix_{sample}`
   * `create_matrix_direct`
     * Step 12 - Create count matrix using original, non-deduplicated reads - this is for showing the "counts" per cell
     * output - `{results}/IsoformMatrix_{sample}_direct`
 * `kristen_isoforms.snake`
   * `make_fastq`
     * Step 14 - Make a fastq files of all reads that have been tagged with the barcode information. Keep the barcode information in the read name
   * `combine_files`
     * Step 14 - Use a custom python script to count up only the reads that mapped uniquely to one isoform
  * `salmon.snake`
    * `sort_minimap_bams`
      * Step 15 - Sort the bam files from `run_minimap_transcriptome` so they can be merged
    * `merge_minimap_bams`
      * Step 15 - Merge the bam files that had been separated for faster processing
    * `split_bams`
      * Step 15 - Split bams based on cell barcode using a custom script. This makes one bamfile per cell
    * `run_salmon`
      * Step 15 - Run salmon on each bam file. This currently runs salmon 9,000+ times, there may be a better way to do this, but I'm not sure yet.
    * `test_salmon`
      * Step 15 - Currently does nothing, but eventually will be an r script that combines the salmon output into one matrix.
  * `full_length_reads.snake`
    * `get_full_length_reads`
      * Step 16 - Makes a bam file for only the full length reads
    * `merge_full_length`
      * Step 16 - Merges the individual bams containing full length reads.

## Input files

## Output files
