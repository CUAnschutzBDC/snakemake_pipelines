# Makes a fastq file with names tagged based on the umi and barcode from scilore
import pysam
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import gzip

infile = snakemake.input[0]
output_file = snakemake.output[0]

samfile = pysam.AlignmentFile(infile, "rb")

def main():
    with gzip.open(output_file, "wt", compresslevel = 6) as out_file:
        for read in samfile:
            gene = "none"

            #pull out tagged info
            for i in read.tags:
                if i[0] == "BC":
                    barcode = i[1]
                if i[0] == "U8":
                    umi = i[1]
                if i[0] == "GE":
                    gene = i[1]
            quality_dict={"phred_quality":read.query_qualities}

            # Reset name to include barcode, umi, and gene
            read_id = read.query_name + "|" + barcode + "|" + umi + "|" + gene
            if read.query_sequence != None:
                # Write record of primary alignements
                record = SeqRecord(
                    Seq(read.query_sequence),
                    id = read_id,
                    name = barcode,
                    description = umi + "_" + gene,
                    letter_annotations = quality_dict,
                    )

                SeqIO.write(record, out_file, "fastq")


if __name__ == '__main__':
    main()