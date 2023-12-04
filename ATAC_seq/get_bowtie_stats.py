import os
import argparse
import glob
import re

def main():
	options = setup()

	if not os.path.exists(options.out_dir):
		os.makedirs(options.out_dir)

	read_files(options.in_dir, options.out_dir, options.trim_method)

#########
# Setup #
#########

def setup():
    """
    Gets command line arguments and returns a Namespace object
    """

    # House keeping to read in arguments from the command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", dest = "in_dir",
        help = "The directory that holds the log files, default is the working directory",
        default = os.getcwd,
        action = "store",
        metavar = "\b")

    parser.add_argument("-o", "--out", dest = "out_dir",
        help = "The directory to put the output, default is outs",
        default = "outs",
        action = "store",
        metavar = "\b")

    parser.add_argument("-t", "--trim-method", dest = "trim_method",
        help = "The trim method used",
        default = "cutadapt",
        action = "store",
        metavar = "\b")

    args = parser.parse_args()

    return(args)

#############
# Functions #
#############

def read_files(in_dir, out_dir, trim_method):
	files = glob.glob(in_dir + "/*.err")
	for file in files:
		sample = re.sub("_" + trim_method + "_trim.err", "", file.split("/")[-1])
		output_file = os.path.join(out_dir, sample + "_bowtie_stats.csv")
		sample_dict = {}
		with open(file, "r") as in_file:
			for line in in_file:
				save_line = line.strip().split(" ")
				if "reads; of these:" in line:
					sample_dict["total_reads"] = save_line[0]
				elif "were paired; of these:" in line:
					sample_dict["paired_reads"] = save_line[0]
					sample_dict["paired_percent"] = save_line[1]
				elif ") aligned concordantly 0 times" in line:
					sample_dict["no_alignment"] = save_line[0]
					sample_dict["no_alignment_percent"] = save_line[1]
				elif "aligned concordantly exactly 1 time" in line:
					sample_dict["one_alignment"] = save_line[0]
					sample_dict["one_alignment_percent"] = save_line[1]
				elif "aligned concordantly >1 times" in line:
					sample_dict["multi_alignment"] = save_line[0]
					sample_dict["multi_alignment_percent"] = save_line[1]
				elif "overall alignment rate" in line:
					sample_dict["overall_alignment"] = save_line[0]
		with open(output_file, "w") as out_file:
			out_file.write("total_reads,paired_reads,paired_percent,no_alignment,no_alignment_percent,one_alignment,one_alignment_percent,multi_alignment,multi_alignment_percent,overall_alignment\n")

			out_file.write("{},{},{},{},{},{},{},{},{},{}\n".format(sample_dict["total_reads"],
				                                                    sample_dict["paired_reads"],
				                                                    sample_dict["paired_percent"],
				                                                    sample_dict["no_alignment"],
				                                                    sample_dict["no_alignment_percent"],
				                                                    sample_dict["one_alignment"],
				                                                    sample_dict["one_alignment_percent"],
				                                                    sample_dict["multi_alignment"],
				                                                    sample_dict["multi_alignment_percent"],
				                                                    sample_dict["overall_alignment"]))


if __name__ == "__main__":
	main()