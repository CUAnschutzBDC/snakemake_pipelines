# Bulk Long Read
A snakemake pipeline that can be used to analyze bulk nanopore long read data. This pipeline aligns reads using minimap2 and then identifies novel transcripts and adds these to a new GTF file and calclulates the abundance of each isoform using [TALON](https://github.com/mortazavilab/TALON).

As of April 2022, the development branch was used because it could take bam files (rather than sam files). This saves both time and space. Running on the released `Talon`, the run on the full promethion dataset continuously failed. The pipeline ran in full using the development branch, but did take over 24 hours.

Benefits of this pipeline are that you can use bulk nanopore data to generate a GTF file with transcripts specific to your cell type, condition, or just not annotated in genocode to use with single cell long read analysis.

TODO: Add the option to quantify with Salmon.
* First, I'll need to create a transcriptome fasta from my Talon created GTF. Can do this with gffread https://www.biostars.org/p/313553/
* Next, I'll need to generate the decoy file with Salmon https://github.com/COMBINE-lab/SalmonTools/blob/master/README.md

* Then I'll need to generate the index and quantify https://combine-lab.github.io/salmon/getting_started/#indexing-txome

Writen by Kristen Wells

## Usage

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install -c conda-forge mamba
mamba install -c bioconda -c conda-forge snakemake
```

3. Update the config file (config.yaml) 
>* RAW_DATA: The path to the directory containing the raw fastq files
>* SAMPLES: The list of sample names. Each name must be unique and must correspond to a directory in the RAW_DATA directory
>* RESULTS: The path to the output directory. It must already exist, you can make it with
```{bash}
mkdir results
```
>* GENOME_FA: The path to the genome fasta file
>* GTF: The path to the GTF associated with the fasta genome
>* GENOME_BUILD: The build of the genome. Will be added to the database files
>* RESULTS: The path to the results directory
>* SICELORE_PATH: The path to the [sicelore package](https://github.com/ucagenomix/sicelore)
>* PICARD_PATH: The path to picard tools
>* MINION_QC_PATH: The path to [minion qc](https://github.com/roblanf/minion_qc)

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`