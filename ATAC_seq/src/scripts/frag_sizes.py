import pysam
import os
import argparse

def main():
	options = setup()

	

	in_dir = "bowtie2_cutadapt_trim"
	in_file = options.sample + "_remove_dup.bam"
	#in_file = options.sample + "_Aligned.sortedByCoord.out.bam"

	out_dir = options.out_dir
	out_file = options.sample + "_frag_sizes.txt"

	in_bam = os.path.join(options.main_dir, in_dir, in_file)
	out_path = os.path.join(options.main_dir, out_dir, out_file)

	samfile = pysam.AlignmentFile(in_bam, "rb")

	read_write_samfile(samfile, out_path)

def setup():
    """
    Gets command line arguments and returns a Namespace object
    """

    # House keeping to read in arguments from the command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", dest = "main_dir",
        help = "The directory that holds the output from the snakefile",
        action = "store",
        metavar = "\b")
    parser.add_argument("-o", "--out-dir", dest = "out_dir",
    	help = "The directory to save results",
    	action = "store",
    	metavar = "\b")
    parser.add_argument("-s", "--sample", dest = "sample",
        help = "the name of the sample",
        default = "WT_mouse1", action = "store",
        metavar = "\b")
   
    args = parser.parse_args()

    return(args)

def read_write_samfile(samfile, out_path):
	with open(out_path, "w") as out_file:
		for read in samfile:
			length = abs(read.template_length)
			out_file.write(str(length) + "\n")


if __name__ == "__main__":
	main()