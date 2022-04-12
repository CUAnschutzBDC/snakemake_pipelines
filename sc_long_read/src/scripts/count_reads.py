# Makes a fastq file with names tagged based on the umi and barcode from scilore
import pysam
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import gzip
import sys

def make_transcript_dict(gtf_file, gene_dict):
    with open(gtf_file, "r") as read_gtf:
        for line in read_gtf:
            if line[0:2] != "##":
                line = line.strip().split("\t")
                if line[2] == "transcript":
                    gene_id = ""
                    transcript_id = ""
                    gene_name = ""
                    info_list = line[8].split(";")
                    for item in info_list:
                        if "gene_id" in item:
                            gene_id = item.split('"')[1]
                        elif "transcript_id" in item:
                            transcript_id = item.split('"')[1]
                        elif "gene_name" in item:
                            gene_name = item.split('"')[1]
                    gene_dict[transcript_id] = gene_name

def read_samfile(samfile, barcode_dict, transcript_dict, results_dict):
    # Make sure I'm correctly reading the last line of the bam
    current_name = ""
    for read in samfile:
        read_name, barcode, umi, gene = read.query_name.split("|")
        barcode_dict[barcode] += 1
        mapping_transcript = read.reference_name
        if mapping_transcript != None:
            gene_name = transcript_dict[mapping_transcript]

            if read.query_name == current_name:
                aligned_list.append(mapping_transcript)
                if gene_name not in current_genes:
                    current_genes.append(gene_name)
            else:
                if current_name != "":
                    if len(current_genes) < 2:
                        if len(aligned_list) == 1:
                            results_dict["|".join(current_genes) + "_" + aligned_list[0]][old_barcode] +=1
                        elif len(aligned_list) > 1:
                            results_dict["|".join(current_genes) + "_undef"][old_barcode] +=1
                        else:
                            sys.exit("No transcripts found")
                old_barcode = barcode
                current_genes = [gene]
                current_name = read.query_name
                aligned_list = [mapping_transcript]
                if gene_name not in gene:
                    current_genes.append(gene_name)

    # Need this for the last line of the file, but should only need it if it passed the if
    # statement above
    if read.query_name == current_name and current_name != "":
        if len(current_genes) < 2:
            if len(aligned_list) == 1:
                results_dict["|".join(current_genes) + "_" + aligned_list[0]][old_barcode] +=1
            elif len(aligned_list) > 1:
                results_dict["|".join(current_genes) + "_undef"][old_barcode] +=1
            else:
                sys.exit("No transcripts found")

def main():
    # Get snakemake info
    input_files = snakemake.input["bam_files"]
    out_file = snakemake.output[0]
    gtf = snakemake.input["gtf"]
    # input_files = "results/minimap2_transcriptome/WT_mouse_1.sorted.bam"
    # out_file = "results/kristen_matrix/WT_mouse_1_counts.txt"
    # gtf = "/beevol/home/rbilab/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf"

    transcript_dict = {}
    make_transcript_dict(gtf, transcript_dict)
    results_dict = defaultdict(lambda : defaultdict(int))

    barcode_dict = defaultdict(int)

    # Determine if the input is a list or string
    if isinstance(input_files, str):
        input_files = [input_files]
    if not isinstance(input_files, list):
        sys.exit("Input files are not in a recognized format!")
        
    for infile in input_files:
        print(infile)
        samfile = pysam.AlignmentFile(infile, "rb")
        read_samfile(samfile, barcode_dict, transcript_dict, results_dict)

    with open(out_file, "w") as output:
        barcode_list = barcode_dict.keys()
        for i in barcode_list:
            output.write("," + i)
        output.write("\n")
        for j in results_dict:
            output.write(j)
            for i in barcode_list:
                output.write("," + str(results_dict[j][i]))
            output.write("\n")


if __name__ == '__main__':
    main()