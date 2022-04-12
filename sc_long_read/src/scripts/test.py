# Makes a fastq file with names tagged based on the umi and barcode from scilore
import pysam
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import gzip

gtf = "/beevol/home/rbilab/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf"


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
    return(gene_dict)
        

def main():
    gene_dict = {}
    transcript_dict = make_transcript_dict(gtf, gene_dict)
    print(gene_dict)


if __name__ == '__main__':
    main()