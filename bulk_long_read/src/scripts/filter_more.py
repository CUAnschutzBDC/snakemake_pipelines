from collections import defaultdict

#filter_file = "filter_transcripts.csv"
filter_file = snakemake.input["filter_file"]

#annotation = "talon_out_talon_read_annot.tsv"
annotation = snakemake.input["annotation"]

#final_filter = "final_filter_transcripts.csv"
final_filter = snakemake.output["filter_file"]

#novel_annotation = "novel_annotation.tsv"
novel_annotation = snakemake.output["annotation"]

filter_dict = defaultdict(int)

write_dict = defaultdict(int)


with open(filter_file, "r") as filtering:
	for line in filtering:
		filter_dict[line.strip()] += 1

with open(annotation, "r") as annotation_file:
	with open(novel_annotation, "w") as novel_write:
		for line in annotation_file:
			if "read_name" in line:
				novel_write.write(line)
			info = line.strip().split("\t")
			gene_ID = info[9]
			transcript_ID = info[10]
			gene_novelty = info[15]
			transcript_novelety = info[16]
			ISM_subtype = info[17]
			dict_key = gene_ID + "," + transcript_ID

			# Keep all known
			if gene_novelty == "Known" and transcript_novelety == "Known":
					write_dict[dict_key] += 1
			if dict_key in filter_dict.keys():
				if gene_novelty == "Known" and transcript_novelety == "NIC":
					write_dict[dict_key] += 1
					if write_dict[dict_key] == 1:
						novel_write.write(line)
				elif gene_novelty == "Known" and transcript_novelety == "NNC":
					write_dict[dict_key] += 1
					if write_dict[dict_key] == 1:
						novel_write.write(line)

with open(final_filter, "w") as write_file:
	for value in write_dict:
		write_file.write(value + "\n")

#print(len(exclude_dict))