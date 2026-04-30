#!/usr/bin/python3
from Bio import SeqIO
import gzip
from io import StringIO
import glob
import os

path_to_files = "/path/to/files/*fastq.gz"

for file in glob.glob(path_to_files):
	filename = os.path.basename(file)

	if filename.endswith(".fastq.gz"):
		prefix = filename[:-9]
		output_file = os.path.join(os.path.dirname(file), prefix + ".fasta")

		with gzip.open(file, "rt") as in_handle, open(output_file, "w") as out_handle:
			count = SeqIO.convert(in_handle, "fastq", out_handle, "fasta")
			print("Converted " + str(count) + " records from " + filename + " to " + prefix + ".fasta")
