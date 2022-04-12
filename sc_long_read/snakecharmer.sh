#! /usr/bin/env bash

#BSUB -J sicelore
#BSUB -o logs/sicelore%J.out
#BSUB -e logs/sicelore%J.err
#BSUB -R "select[mem>4] rusage[mem=4]" 
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

module load java/1.8
module load samtools
#module load ucsc/v308
module load R/4.0.3
 
# Function to run snakemake
run_snakemake() {
    local num_jobs=$1
    local config_file=$2

    args='
        -o {log}.out 
        -e {log}.err 
        -J {params.job_name} 
        -R "{params.memory} span[hosts=1] " 
        -n {threads} 
        -q rna '

    snakemake \
        --snakefile Snakefile \
        --drmaa "$args" \
        --jobs $num_jobs \
        --latency-wait 60 \
        --rerun-incomplete \
        --configfile $config_file
}

run_snakemake 24 config.yaml
