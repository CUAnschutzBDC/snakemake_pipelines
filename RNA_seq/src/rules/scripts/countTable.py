from collections import defaultdict
import re

file_list = snakemake.input
output_file = snakemake.output[0]

trim_method = snakemake.params.trim_method
results = snakemake.params.results

gene_dict = defaultdict(list)
sample_list = list()

file_path = results + "/featureCount" + "_" + trim_method + "_trim/"

for i in file_list:
	sample_name = re.sub(rf'{file_path}', '', i)
	sample_name = re.sub(r'_countsOutput', '', sample_name)
	sample_list.append(sample_name)
	with open(i, "r") as countFile:
		for line in countFile:
			line = line.strip().split("\t")
			if "ENSMUSG" in line[0]:
				ens_id = line[0]
				gene = line[6]
				count = line[8]
				gene_ens = gene + "_" + ens_id
				gene_dict[gene_ens].append(count)

with open(output_file, "w") as count_file:
	count_file.write("gene" + "\t" + "\t".join(sample_list) + "\n")
	for gene in gene_dict:
		count_file.write(gene + "\t" + "\t".join(gene_dict[gene]) + "\n")