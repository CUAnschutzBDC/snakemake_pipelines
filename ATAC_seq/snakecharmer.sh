#!/bin/bash 
#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=1gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=christopher.2.hill@cuanschutz.edu

set -o nounset -o pipefail -o errexit -x
mkdir -p logs
module load anaconda
conda activate snakemake8
#export APPTAINER_BIND="/scratch/alpine:/scratch/alpine,/pl/active/Anschutz_BDC:/pl/active/Anschutz_BDC"

# Run snakemake pipeline
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 60 \
    --latency-wait 60 \
    --rerun-incomplete \
    --workflow-profile profiles/default \
    --ignore-incomplete
    