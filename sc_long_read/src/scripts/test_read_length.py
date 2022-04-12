import pysam
from collections import defaultdict

length_dict_all = defaultdict(lambda: defaultdict(int))
length_dict_umi = defaultdict(lambda: defaultdict(int))

# Read in files from snakemake
all_file = snakemake.input["all_input"]
umi_file = snakemake.input["umi_input"]
outfile_all_gene = snakemake.output["all_output_gene"]
outfile_umi_gene = snakemake.output["umi_output_gene"]
outfile_all_no_gene = snakemake.output["all_output_no_gene"]
outfile_umi_no_gene = snakemake.output["umi_output_no_gene"]

#umi_file = ["0016.GEUS10xAttributes.sorted_umifound_.bam"]

#all_file = ["0016.GEUS10xAttributes.sorted.bam"]

# Find gene lengths from all reads
for file_name in umi_file:
    samfile = pysam.AlignmentFile(file_name, "rb")
    for read in samfile.fetch():
        sequence_len = len(read.get_tag("US"))
        if ('BG', '') in read.get_tags():
            length_dict_umi["no_gene"][sequence_len] += 1
        else:
            length_dict_umi["gene"][sequence_len] += 1



# Find gene lengths from kept reads
for file_name in all_file:
    samfile = pysam.AlignmentFile(file_name, "rb")
    for read in samfile.fetch():
        sequence_len = len(read.get_tag("US"))
        if ('BG', '') in read.get_tags():
            length_dict_all["no_gene"][sequence_len] += 1
        else:
            length_dict_all["gene"][sequence_len] += 1


with open(outfile_all_gene, "w") as outfile:
    outfile.write("{},{}\n".format("length", "count"))
    for i in length_dict_all["gene"]:
        outfile.write("{},{}\n".format(i, length_dict_all["gene"][i]))

with open(outfile_all_no_gene, "w") as outfile:
    outfile.write("{},{}\n".format("length", "count"))
    for i in length_dict_all["no_gene"]:
        outfile.write("{},{}\n".format(i, length_dict_all["no_gene"][i]))

with open(outfile_umi_gene, "w") as outfile:
    outfile.write("{},{}\n".format("length", "count"))
    for i in length_dict_umi["gene"]:
        outfile.write("{},{}\n".format(i, length_dict_umi["gene"][i]))

with open(outfile_umi_no_gene, "w") as outfile:
    outfile.write("{},{}\n".format("length", "count"))
    for i in length_dict_umi["no_gene"]:
        outfile.write("{},{}\n".format(i, length_dict_umi["no_gene"][i]))