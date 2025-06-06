#!/bin/bash 

#SBATCH --job-name=sif
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=6gb
#SBATCH --output=logs/sif.out
#SBATCH --partition=amilan

apptainer build atac_chip_v1.sif Singularity