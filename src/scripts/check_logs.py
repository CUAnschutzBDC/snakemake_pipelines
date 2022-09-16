import re
import sys
import argparse
import glob
import os


def main():
	options = setup()
	failed_list = []

	if options.in_dir == "none":
		print("no directory provided, looking in current directory")
	else:
		os.chdir(options.in_dir)

	check_files(failed_list)
	print_failures(failed_list)

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
        help = "The directory that holds the log files",
        default = "none",
        action = "store",
        metavar = "\b")

    args = parser.parse_args()

    return(args)

def check_files(failed_list):
	"""
	Goes through all out log files and looks for the "successfully completed."
	message. If this isn't found, the file name is added to the failed list.
	"""
	for file in glob.glob("*.out"):
		success = False
		with open(file, "r") as input_file:
			for line in input_file:
				if re.search("Successfully completed.", line):
					success = True
					break
		if not success:
			failed_list.append(file)

def print_failures(failed_list):
	"""
	Prints out all files to the terminal.
	"""
	if len(failed_list) == 0:
		print("No rules failed!")
	else:
		print("Failed files:")
		for failed_file in failed_list:
			print(failed_file)

if __name__ == "__main__":
	main()