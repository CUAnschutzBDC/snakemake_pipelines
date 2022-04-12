import pysam

samfile = pysam.AlignmentFile(infile, "rb")

with pysam.AlignmentFile(output_file, "wb", template = samfile) as out_file:
	for read in samfile:
		if abs(read.reference_length) < 300:
			out_file.write(read)