# ChIP_seq analysis pipeline

A snakemake pipeline that can be used to run bulk ChIP-seq analysis. Trimming can be done with bbduk, cutadapt, or no triming. Outputs fastqc summary files and peak files from Macs2 that can be analyzed using the rmd script

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
>* SAMPLE_TABLE: path to a file consisting of at least five columns including all sampels and controls.
>   * The first column should be titled `sample` and contain the name of the sample.
>   * The second should be titled `fastq1` and contain the path to the associated fastq file
>   * The third *optional* column should be titled `fastq2` and be used if you have paired end data. This should contain the path to the associated read 2 fastq file.
>   * The fourth should be titled data_dir and contain the path to the data
>   * The fifth should be labeled control and contain name of the associated control file (or "None" if the file was a control)
>   * The sixth should be labeled replicate and contain the replicate number of the sample this is for use by ChIPQC
>   * The seventh *optional* labeled Tissue can contain the name of the tissue. This is for use by ChIPQC
>   * The eigth *optional* labeled factor can list the TF names. This is for use by ChIPQC
>   * The ninth *optional* labeled contidion can have information about your sample. This is for use by ChIPQC.
>* PROJECT: The name of the project, will be the name of the output counts matrix
>* GENOME: The path to the genome directory created by star
>* Annotation: The annotation used "mm10" for example
>* RESULTS: The path to the results directory
>* ADAPTORS: *Optional*, only include if using bbduk for adaptor trimming. Path to the adaptors file in the bbtools package
>* TRIM_METHOD: What trim method to use. Can be "bbduk", "cutadapt", or "no"
>* GENOME_SIZE: The size of the genome. Some examples (taken from MACS2) are in the configfile
>* SANDBOX: Path to sandbox to upload data.
>* PE and SE: extra paramamters for all of the jobs run

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`