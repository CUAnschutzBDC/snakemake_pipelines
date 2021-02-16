# Analysis pipelines
Writen by Kristen Wells

A collection of snakemake pipelines to analyze RNA datasets. Current pipelines include:

* A scRNA-sequencing pipeline that analyzed 10x genomics datasets. Can be used with scRNA-seq, scCITE-seq, scVDJ-seq, and hashtagging (and any combination of those).
* A basic bulk RNA-sequencing pipeline that includes options for trimming (cutadapt, bbduk or none). The pipeline includes alignment with star and read counting with FeatureCounts.
* A pipieline to subset fastq files into smaller chunks.

Pipelines that will be added soon include:
* A bulk RNA-seq pipeline that is capable of analyzing NET seq data with or without spikie ins.


All snakemake pipelines presented here work best when run on an lsf cluster.

Snakemake must be installed to use these pipelines.

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install snakemake -c bioconda -c conda-forge
```

Steps to run indivdual pipelines are within their own directories.