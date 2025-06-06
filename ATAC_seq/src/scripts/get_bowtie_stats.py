import os
import argparse
import glob
import re

def main():
	options = setup()

	if not os.path.exists(options.out_dir):
		os.makedirs(options.out_dir)

	out_file = os.path.join(options.out_dir, "bowtie2_stats.csv")
	full_dict = read_files(options.in_dir, options.trim_method)
	write_output(full_dict, out_file)

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

def read_files(in_dir, trim_method):
	"""
	This reads the log files output by bowtie2 to determine the alignment rates of reads
	
	The output is a layered dictionary. The inner layer is information per sample,
	all of these inner dicts are layered into a full dict where sample is the key.
	"""
	files = glob.glob(in_dir + "/*.err")
	full_dict = {}
	for file in files:
		sample = re.sub("_" + trim_method + "_trim.err", "", file.split("/")[-1])
		sample = re.sub("bowtie2_", "", sample)
		#output_file = os.path.join(out_dir, sample + "_bowtie_stats.csv")
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
		full_dict[sample] = sample_dict
	return(full_dict)

def write_output(full_dict, output_file):
	"""
	Here the dictionary created by read_files is written to an output file.
	"""
	with open(output_file, "w") as out_file:
		out_file.write("sample,total_reads,paired_reads,paired_percent,no_alignment,no_alignment_percent,one_alignment,one_alignment_percent,multi_alignment,multi_alignment_percent,overall_alignment\n")

		for sample in full_dict:
			sample_dict = full_dict[sample]
			out_file.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(
				sample, sample_dict["total_reads"], sample_dict["paired_reads"],
				sample_dict["paired_percent"], sample_dict["no_alignment"],
				sample_dict["no_alignment_percent"], sample_dict["one_alignment"],
				sample_dict["one_alignment_percent"], sample_dict["multi_alignment"],
				sample_dict["multi_alignment_percent"], sample_dict["overall_alignment"]))


if __name__ == "__main__":
	main()