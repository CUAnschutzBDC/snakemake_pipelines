"""
Count tags from the scilore pipeline in all reads and in reads thrown out to
get an idea for why reads are thrown out.
"""

import pysam
from collections import defaultdict

barcode_dict_all = defaultdict(int)
barcode_dict_umi = defaultdict(int)

# Read in files from snakemake
all_file = snakemake.input["all_input"]
umi_file = snakemake.input["umi_input"]
outfile = snakemake.output[0]

# Set reads to 0
reads_all = 0
reads_umi = 0

# Find percentages from all reads
for file_name in umi_file:
    samfile = pysam.AlignmentFile(file_name, "rb")
    for read in samfile.fetch():
        all_tags = read.get_tags()
        for tag in all_tags:
            barcode_dict_umi[tag[0]] += 1
        reads_umi += 1


# Find percentages from kept reads
for file_name in all_file:
    samfile = pysam.AlignmentFile(file_name, "rb")
    for read in samfile.fetch():
        all_tags = read.get_tags()
        for tag in all_tags:
            barcode_dict_all[tag[0]] += 1
        reads_all += 1

with open(outfile, "w") as out:
    out.write("{}\t{}\t{}\n".format("tag", "count", "percent"))
    for tag in barcode_dict_umi:
        # Calculate the percent of reads in each category
        percent = barcode_dict_umi[tag]/reads_umi
        out.write("{}\t{}\t{}\n".format(tag, barcode_dict_umi[tag], percent))
    out.write("\n")
    out.write("{}\t{}\t{}\n".format("tag", "count", "percent"))
    for tag in barcode_dict_all:
        percent = barcode_dict_all[tag]/reads_all
        out.write("{}\t{}\t{}\n".format(tag, barcode_dict_all[tag], percent))