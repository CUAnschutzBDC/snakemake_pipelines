# subset_fastq
A snakemake pipeline that can be used to subset fastq files. This pipeline works best when run on an lsf cluster.

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
>* data_dir: Location of the raw data 
>* samples: the list of samples you want to test 

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`


