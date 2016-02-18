#############################
#18th Feb 2016				#
#dr.mark.schultzm@gmail.com #
#github: schultzm			#
#############################

"""
Scans each line of an alleles.csv file (csv, header == strain names,
column 1 == snp pos, subsequent columns are allele calls) and checks if snp pos
is between ranges in regions file (no header, tab delimited). If it is,
outputs the line of the alleles.csv file to subset file.  Always outputs
header (strain names) of Alleles.csv file.
"""
import argparse
import os

#set up the arguments parser to deal with the command line input
parser = argparse.ArgumentParser(description = "Subset an Alleles.csv file by \
									pulling out regions between start and stop\
									coordinates in tab delimited text file (no\
									header).")
parser.add_argument('-i', '--input', help = 'Input alleles.csv to pull the \
					sub-alleles.csv from', required=True)
parser.add_argument('-r', '--regions', help = 'Regions to target', \
					required = True)
ARGS = parser.parse_args()

def get_regions(alleles_csv, regions_file):
	"""
	Main function
	"""
	input_handle_name = os.path.basename(alleles_csv)
	input_handle_prefix = os.path.splitext(input_handle_name)[0]
	regions = open(regions_file, "r")
	blocks = []
	for line in regions:
		positions = [i for i in line.rstrip().split('\t')]
		positions = map(int, positions)
		blocks.append(positions)
	output_handle = open(input_handle_prefix+'_SNPsitesSubset.csv', "w")
	with open(alleles_csv, "r") as input_handle:
		output_handle.write(next(input_handle))
		for line in input_handle:
			snp_pos = int(line.rstrip().split(',')[0])
			#check if snp_pos between range of internal lists in blocks list
			if len(filter(lambda x: x[0] <= snp_pos <= x[1], blocks)) > 0:
				output_handle.write(line)

get_regions(ARGS.input, ARGS.regions)
