# scRNA_seq
A snakemake pipeline that can be used to analyze the scRNA-seq data generated for the publication [Calcium-dependent transcriptional changes in human pancreatic islet cells reveal functional diversity in islet cell subtypes](https://link.springer.com/article/10.1007/s00125-022-05718-1)

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
>* RAW_DATA: Location of the raw data 
>* SAMPLES: the list of samples you want to test. This is the name that will be in the output files. The order must be identical to the order of RNA_SAMPLES, ADT_SAMPLES, and VDJ_SAMPLES
>* AGGR_SAMPLES: Which samples should be aggregated together at the end of the pipeline. These should be samples that were split between multiple 10x runs.
>* RNA_SAMPLES: The name of the RNA samples. Not the full fastq name, but the name that is the same for all samples. Likely ends in "GEX"
>* ADT_SAMPLES: *optional* The name of the samples from CITE-seq or hashtagging. If CITE-seq and hashtagging files are separate, include them both separated by a comma. Put samples in the same order as their RNA counterparts. If CITE-seq or hashtagging were not performed, leave this blank.
>* VDJ_SAMPLES: *optional* The name of the samples from VDJ-seq. Put samples in the same order as their RNA counterparts. If VDJ-sequencing was not peformed, leave this blank.
>* RESULTS: Path to the output directory
>* GENOME: Path to the cellranger genome reference
>* ADT_REF: *optional* Path to the ADT-reference. Should be a comma separated file with the columns described in the 10x tutorial: id, name, read, pattern, sequence, and feature_type. The feature_type will be Antibody Capture. The name will be the name in the final output matrix. Leave this blank if CITE-seq or hashtagging were not performed.
>* VDJ_REF: Path to the cellranger VDJ reference. If VDJ sequencing were not performed, leave this blank.
>* MAX_JOBS: The maximum number of jobs that can be submitted by cell ranger at a time
>* LSF_TEMPLATE: Path to an LSF template. One is included in this git repo.
>* CHEMISTRY: *optional* Arguments to the `--chemstiry` flag in cellranger count. If left blank, chemistry will be `auto`. Only use this if the pipeline failed because of a failure to detect the chemistry. Can be filled in for only some samples.

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`

6. I highly recommend looking at the csv files that are generated and passed to cell ranger to ensure that the correct fastq files have been detected for each sample.

7. Follow the analysis scripts in the order shown in the `src` folder. Analysis requires the `scAnalysisR` package available on `github`:

```R
install.packages("devtools")
library(devtools)
install_github("CUAnschutzBDC/scAnalysisR")
```
