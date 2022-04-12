import pysam
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import gzip

#samfile = pysam.AlignmentFile("0001.GEUS10xAttributes.sorted_umifound_.bam", "rb")

input_files = snakemake.input["bam_files"]
results_path = snakemake.params["results_dir"]


def read_samfile(samfile, output_dict):
    for read in samfile:
        gene = "none"
        for i in read.tags:
            if i[0] == "BC":
                barcode = i[1]
            if i[0] == "U8":
                umi = i[1]
            if i[0] == "IG":
                gene = i[1]
        quality_dict={"phred_quality":read.query_qualities}
        file_name = results_path + barcode + ".fastq"
        #print(read.qual)
        if read.query_sequence != None:
            record = SeqRecord(
                Seq(read.query_sequence),
                id = read.query_name,
                name = barcode,
                description = umi + "_" + gene,
                letter_annotations = quality_dict,
                )

            output_dict[file_name].append(record)

def main():
    for infile in input_files:
        # initialize output dictionary
        output_dict = defaultdict(list)
        print(infile)

        samfile = pysam.AlignmentFile(infile, "rb")
        read_samfile(samfile, output_dict)
    
        for file in output_dict:
            print(file)
            with open(file, "a") as out_file:
                for sequence in output_dict[file]:
                    SeqIO.write(sequence, out_file, "fastq")


if __name__ == '__main__':
    main()