# ATAC-seq analysis pipeline

A snakemake pipeline that can be used to run bulk ATAC-seq analysis. Can chose between cutadapt, bbduk or no adapter trimming. Outputs fastqc summary files, bowtie2 summary files, and a counts matrix that can be analyzed using the rmd script

Writen by Kristen Wells

To use:

1. Download and install `miniconda3`: For Linux
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install `Snakemake` and make sure `singularity` is available on your system:
```bash
conda install snakemake -c bioconda -c conda-forge
```

3. Update the config file (config.yaml) 
>* SAMPLE_TABLE: path to a file consisting of at least two columns.
>   * The first column should be titled `sample` and contain the name of the sample.
>   * The second should be titled `fastq1` and contain the path to the associated fastq file
>   * The third *optional* column should be titled `fastq2` and be used if you have paired end data. This should contain the path to the associated read 2 fastq file.
>   * The fourth *optional* column is named `spike-in` with TRUE/FALSE values per sample indicating if spikeins were used. If this column is not included it will be assumed that spike-ins were not included.
>* PROJECT: The name of the project, will be the name of the output counts matrix
>* GENOME: The path to the genome directory created by bowtie2
>* GTF: The path to the GTF associated with the star genome
>* RESULTS: The path to the results directory
>* ADAPTORS: *Optional*, only include if using bbduk for adaptor trimming. Path to the adaptors file in the bbtools package
>* TRIM_METHOD: What trim method to use. Can be "bbduk", "cutadapt", or "no"
>* PE and SE: extra paramamters for all of the jobs run
>* BLACKLIST: Path to the blacklist file. Can download [here](https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/)
>* HMMRATAC_FILE: The path to the jar file for hmmratac. Can download [here](https://github.com/LiuLabUB/HMMRATAC/releases)
>* GENERAL_CONTAINER: The path to the general container. Can be built using the recipe in `docker/general_docker` (or eventually downloaded)
>* PICARD_CONTAINER: The path to the picard container. Can be built using the recipe in `docker/picard_docker` (or eventually downloaded)
>* PICARD_JAR: Path to the picard jar. This is included in the repo so the path should not need to be changed.
>* R_CONTINER: Path to the r docker container. Can be built using the recipe in `docker/r_docker` (or eventually downloaded).
>* SEQ_INFO: Information about your species from bioconductor. This is because the txdb package often failed for me. Made by running `seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)` in R.
>* TRANSCRIPTS: Path to transcripts from bioconductor. This is because the txdb package often failed for me. Made by running `txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)`
>* BSGENOME: Path to the bsgenome. Can be generated by running `genome <- BSgenome.Musmusculus.UCSC.mm10::Mmusculus`



4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`


## Bowtie2 arguments

```
--very-sensitive
Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50



-N <int>
Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity. Default: 0.

-L <int>
Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 22 and 20 in --end-to-end mode and in --local mode.

-i <func>
Sets a function governing the interval between seed substrings to use during multiseed alignment. For instance, if the read has 30 characters, and seed length is 10, and the seed interval is 6, the seeds extracted will be:

Read:      TAGCTACGCTCTACGCTATCATGCATAAAC
Seed 1 fw: TAGCTACGCT
Seed 1 rc: AGCGTAGCTA
Seed 2 fw:       CGCTCTACGC
Seed 2 rc:       GCGTAGAGCG
Seed 3 fw:             ACGCTATCAT
Seed 3 rc:             ATGATAGCGT
Seed 4 fw:                   TCATGCATAA
Seed 4 rc:                   TTATGCATGA
Since it's best to use longer intervals for longer reads, this parameter sets the interval as a function of the read length, rather than a single one-size-fits-all number. For instance, specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. See also: setting function options. If the function returns a result less than 1, it is rounded up to 1. Default: the --sensitive preset is used by default, which sets -i to S,1,1.15 in --end-to-end mode to -i S,1,0.75 in --local mode.

-D <int>
Up to <int> consecutive seed extension attempts can "fail" before Bowtie 2 moves on, using the alignments found so far. A seed extension "fails" if it does not yield a new best or a new second-best alignment. This limit is automatically adjusted up when -k or -a are specified. Default: 15.

-R <int>
<int> is the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds. When "re-seeding," Bowtie 2 simply chooses a new set of reads (same length, same number of mismatches allowed) at different offsets and searches for more alignments. A read is considered to have repetitive seeds if the total number of seed hits divided by the number of seeds that aligned at least once is greater than 300. Default: 2.



-k <int>
By default, bowtie2 searches for distinct, valid alignments for each read. When it finds a valid alignment, it continues looking for alignments that are nearly as good or better. The best alignment found is reported (randomly selected from among best if tied). Information about the best alignments is used to estimate mapping quality and to set SAM optional fields, such as AS:i and XS:i.

When -k is specified, however, bowtie2 behaves differently. Instead, it searches for at most <int> distinct, valid alignments for each read. The search terminates when it can't find more distinct valid alignments, or when it finds <int>, whichever happens first. All alignments found are reported in descending order by alignment score. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM 'secondary' bit (which equals 256) set in its FLAGS field. For reads that have more than <int> distinct, valid alignments, bowtie2 does not guarantee that the <int> alignments reported are the best possible in terms of alignment score. -k is mutually exclusive with -a.

Note: Bowtie 2 is not designed with large values for -k in mind, and when aligning reads to long, repetitive genomes large -k can be very, very slow.


-X/--maxins <int>
The maximum fragment length for valid paired-end alignments. E.g. if -X 100 is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). A 61-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the -X constraint is applied with respect to the untrimmed mates, not the trimmed mates.

The larger the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences between -I and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient.

Default: 500.


--rg-id <text>
Set the read group ID to <text>. This causes the SAM @RG header line to be printed, with <text> as the value associated with the ID: tag. It also causes the RG:Z: extra field to be attached to each SAM output record, with value set to <text>.

--rg <text>
Add <text> (usually of the form TAG:VAL, e.g. SM:Pool1) as a field on the @RG header line. Note: in order for the @RG line to appear, --rg-id must also be specified. This is because the ID tag is required by the SAM Spec. Specify --rg multiple times to set multiple fields. See the SAM Spec for details about what fields are legal.
```