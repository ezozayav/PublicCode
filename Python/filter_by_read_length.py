##########################
#dr.mark.schultz@gmail.com
#github schultzm
#date 06/11/15
##########################

import os
from Bio import SeqIO
import argparse

#set up the argument parser
parser = argparse.ArgumentParser(description = "Will filter seqs in file to be above a certain length.")
parser.add_argument('-i', '--input_contigs', nargs ='+', help = 'Input copntig (seq) files (.fasta or mfasta) to filter', required=True)
parser.add_argument('-m', '--min_length', type=int, help = 'Minimum length of contig to pass filter (integer)', default = 100)
args = parser.parse_args()

#This function does all the work
def filter_contigs(infile, length):
	"""
	Input file to filter contigs from
	"""
	short_sequences = [] # Setup an empty list
	with open(infile, "r") as input_handle:
		with open(infile.replace('mfasta', '')+'greaterthan'+str(length)+'bp.mfasta', 'w') as output_handle:
			for record in SeqIO.parse(input_handle, "fasta"):
				if len(record.seq) > length :
					# Add this record to our list
					SeqIO.write(record, output_handle, "fasta")

#execute the filter_contigs function for each contig file
for i in args.input_contigs:
	filter_contigs(i, args.min_length)


