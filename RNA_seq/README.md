# RNA_seq analysis pipeline

A snakemake pipeline that can be used to run bulk RNA-seq analysis. Can chose between cutadapt, bbduk or no adapter trimming. Outputs fastqc summary files, star summary files, and a counts matrix that can be analyzed using the rmd script

Writen by Kristen Wells

To use:

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install snakemake -c bioconda -c conda-forge
```

3. Update the config file (config.yaml) 
>* SAMPLE_TABLE: path to a file consisting of at least two columns.
>   * The first column should be titled `sample` and contain the name of the sample.
>   * The second should be titled `fastq1` and contain the path to the associated fastq file
>   * The third *optional* column should be titled `fastq2` and be used if you have paired end data. This should contain the path to the associated read 2 fastq file.
>   * The fourth *optional* column is named `spike-in` with TRUE/FALSE values per sample indicating if spikeins were used. If this column is not included it will be assumed that spike-ins were not included.
>* PROJECT: The name of the project, will be the name of the output counts matrix
>* GENOME: The path to the genome directory created by star
>* COMBINED_GENOME: *Optional*, only include this if spike-ins were used. This is the genome of both the main organism and spike-in organism
>* GTF: The path to the GTF associated with the star genome
>* COMBINED_GTF: *Optional*, only include with spike-ins, the path to the combined GTF file
>* RESULTS: The path to the results directory
>* ADAPTORS: *Optional*, only include if using bbduk for adaptor trimming. Path to the adaptors file in the bbtools package
>* TRIM_METHOD: What trim method to use. Can be "bbduk", "cutadapt", or "no"
>* PE and SE: extra paramamters for all of the jobs run

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`