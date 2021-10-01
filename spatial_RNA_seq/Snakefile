""" Snake pipeline for running spaceranger """

# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")
#import subprocess
import glob
import os 
#import re
#from collections import defaultdict

# Parameters from config.yaml
RAW_DATA      = config["RAW_DATA"]
SAMPLES       = config["SAMPLES"]
RESULTS       = config["RESULTS"]
TRANSCRIPTOME = config["TRANSCRIPTOME"]
IMAGES        = config["IMAGES"]
LSF_TEMPLATE  = config["LSF_TEMPLATE"]
PROBE_SET     = config["PROBE_SET"]
SLIDE_PATH    = config["SLIDE_PATH"]


# Function to check paths for input files/directories
def _check_path(path):
    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")


if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)
    
if PROBE_SET:
    FFPE = True
else:
    FFPE = False

# Pull out sample names
sample_names = SAMPLES.keys()

# set fastq information
FASTQ_INFO      = "_S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz"
FASTQ_DIR = RESULTS + "/fastqs"
if not os.path.exists(FASTQ_DIR):
    os.makedirs(FASTQ_DIR)

# Check dirs
RESULTS  = _check_path(RESULTS)
TRANSCRIPTOME = _check_path(TRANSCRIPTOME)
if FFPE:
    PROBE_SET = _check_path(PROBE_SET)

if LSF_TEMPLATE:
    LSF_TEMPLATE = _check_path(LSF_TEMPLATE)
else:
    LSF_TEMPLATE = "lsf"

SLIDE_PATH = _check_path(SLIDE_PATH)

rule all:
    input:
        # Merge fastq files
        expand(
        	"{results}/logs/merge_fastqs/{sample}_merge_fastqs_done.out",
        	results = RESULTS, sample = sample_names
        	),
        expand(
            "{results}/logs/run_spaceranger/{sample}_spaceranger_done.out",
            results = RESULTS, sample = sample_names
            )

include: "src/rules/spaceranger.snake"