# This script takes the bam file aligned by minimap 2 and contains UMI and barcode information
# in the read header (put there by make_fastq.py) and separates the reads into individual
# bam files for each cell.

import pysam
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import gzip

#samfile = pysam.AlignmentFile("0001.GEUS10xAttributes.sorted_umifound_.bam", "rb")

infile = snakemake.input["bam_file"]
results_path = snakemake.params["out_dir"]

if not os.path.exists(results_path):
    os.makedirs(results_path)

#infile = "/beevol/home/wellskri/Analysis/Lori_Sussel/Maria_Hansen/long_read_sc/nanopore_new_gtf/results/minimap2_transcriptome/WT_mouse_1/0001.sorted.bam"
#results_path = "/beevol/home/wellskri/Analysis/Lori_Sussel/Maria_Hansen/long_read_sc/nanopore_new_gtf/test"

def read_samfile(samfile, output_dict):
    for read in samfile:
        gene = "none"
        read_name, barcode, umi, gene = read.query_name.split("|")
        file_name = os.path.join(results_path, barcode + ".bam")
        output_dict[file_name].append(read)

def main():
    # initialize output dictionary
    output_dict = defaultdict(list)
    print(infile)

    samfile = pysam.AlignmentFile(infile, "rb")
    read_samfile(samfile, output_dict)

    for file in output_dict:
        with pysam.AlignmentFile(file, "wb", template = samfile) as out_file:
            for sequence in output_dict[file]:
                out_file.write(sequence)


if __name__ == '__main__':
    main()