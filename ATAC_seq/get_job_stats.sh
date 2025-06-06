#!/bin/bash 

log_dir=logs
script="/pl/active/Anschutz_BDC/resources/snakemake/snakemake_pipelines/src/scripts/job_efficiency.py"

# All jobs
python $script \
	-d $log_dir \
	-o job_stats_all.csv \
	--keep_failed

# Only completed jobs
python $script \
	-d $log_dir \
	-o job_stats.csv