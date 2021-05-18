from collections import defaultdict
import re

file_list = snakemake.input
output_file = snakemake.output[0]

trim_method = snakemake.wildcards.trim_method
results = snakemake.wildcards.results

stat_dict = defaultdict(list)
sample_list = list()

file_path = results + "/star" + "_" + trim_method + "_trim/"

for i in file_list:
	sample_name = re.sub(rf'{file_path}', '', i)
	sample_name = re.sub(r'_Log\.final\.out', '', sample_name)
	sample_list.append(sample_name)
	with open(i, "r") as countFile:
		for line in countFile:
			print(line)
			if "READS" not in line:
				line = line.strip().split(" |\t")
				if len(line) == 2:
					stat_type, value = line
					stat_dict[stat_type].append(value)

with open(output_file, "w") as stat_file:
	stat_file.write("statistic" + "\t" + "\t".join(sample_list) + "\n")
	for statistic in stat_dict:
		stat_file.write(statistic + "\t" + "\t".join(stat_dict[statistic]) + "\n")