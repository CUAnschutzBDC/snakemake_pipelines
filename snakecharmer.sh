#!/usr/bin/env bash

#BSUB -J RNAseq
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

# Load modules
# Use locally installed packages for Majiq
# module load fastqc/0.11.7
# module load samtools/1.12
# module load STAR/2.7.9a
# module load subread/1.6.2
# module load bbtools/38.86
module load singularity/3.9.2

# LSF arguments
args=' 
  -q rna 
  -o {log}.out 
  -e {log}.err 
  -J {params.job_name} 
  -R "{params.memory} span[hosts=1] " 
  -n {threads} ' 

# Run snakemake pipeline
snakemake \
    --drmaa "$args" \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 60 \
    --latency-wait 60 \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args "--bind /beevol/home/rbilab --bind /beevol/home/wellskri --bind /beevol/home/wellskri/packages/bin