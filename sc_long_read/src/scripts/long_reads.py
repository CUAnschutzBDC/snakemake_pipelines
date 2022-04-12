# This script takes the bam file aligned by minimap 2 and contains UMI and barcode information
# in the read header (put there by make_fastq.py) and keeps only reads that map to a certain 

# TODO - write so it will take in a list of files

import pysam
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import gzip

# infile = snakemake.input[0]
# output_file = snakemake.output[0]
# percent_cutoff = snakemake.params["cutoff"]

infile = "results/tag_UMI_consensus/WT_mouse_1.GEBC.sorted.bam"
output_file = "results/gviz/WT_mouse_1_full_length.bam"
percent_cutoff = 0.8

def write_long_reads(samfile, output_file, length_dict):
    with pysam.AlignmentFile(output_file, "wb", template = samfile) as out_file:
        current_name = ""
        sequence = ""
        for read in samfile:
            transcript = read.reference_name
            if transcript != None:
                transcipt_length = length_dict[transcript]
                # Pull out sequence information if it's the first time seeing the read
                if read.query_name != current_name:
                    current_name = read.query_name
                    sequence = Seq(read.query_sequence)
                # Find the percent of the full transcript length
                sequence_percent = len(sequence)/transcipt_length
                if sequence_percent > percent_cutoff:
                    out_file.write(read)


def make_length_dict(samfile, length_dict):
    length_list = samfile.header["SQ"]
    for i in length_list:
        length_dict[i["SN"]] = i["LN"]


def main():

    samfile = pysam.AlignmentFile(infile, "rb")
    length_dict = {}
    make_length_dict(samfile, length_dict)
    write_long_reads(samfile, output_file, length_dict)



if __name__ == '__main__':
    main()