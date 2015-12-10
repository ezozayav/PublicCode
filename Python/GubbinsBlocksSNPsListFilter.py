#Mark Schultz 041214
#dr.mark.schultz@gmail.com

import argparse
import sys

parser = argparse.ArgumentParser(description = "Provide SNP list and Gubbins output blocks file.  For each SNP in list, will scan blocks file to see if SNP falls within a recombinant block.  If True, the SNP will be excluded from the output list of SNPs (or included if 'recombinant' SNPs are requested, specified on '-f' switch).")
parser.add_argument('-s', '--snp_list', help = "File containing SNP genome coordinates (rows) in a single column. No header row.", required = True)
parser.add_argument('-b', '--blocks_file', help = "Blocks file output from Gubbins, with strain name in column 1, start block position in column 2 and end block in column 3. Tab-delimited, no header line in file.", required=True)
parser.add_argument('-f', '--filter_option', help = "Enter 'recombinant' or 'nonrecombinant' depending on whether you want the output file to include the recombinant (SNPs that fall within the Gubbins blocks) or the non-recombinant SNPs (those that don't fall within the Gubbins blocks) in the input_list. Default = 'nonrecombinant'.", default = 'nonrecombinant', required=False)
parser.add_argument('-p', '--output_filename_prefix', help = "Some meaningful text to append to the start of the output snp_list.", required=True)
args = parser.parse_args()

print ""

#get the unique occurences of block ranges in the blocks file
def get_blocks(blocks_file):
	blocks = []
	blocks_unique = []
	with open(blocks_file, 'r') as blocks_handle:
		for line in blocks_handle:
			x = line.rstrip("\n").replace(",","\t").split("\t")
			blocks.append([int(x[1]),int(x[2])])
	blocks.sort()
	for i in blocks:
		if i not in blocks_unique:
			blocks_unique.append(i)
	blocks_unique.sort
	return blocks_unique

def get_snps(snp_list):
	snps = []
	with open(snp_list, 'r') as snps_handle:
		for line in snps_handle:
			y = line.rstrip("\n")
			snps.append(int(y))
	snps.sort()
	return snps


def filter_snps(blocks, snps, filter, prefix):
	snps_in_blocks = []
	for i in blocks:
		for j in snps:
			if j >= i[0] and j <= i[1]:
				if j not in snps_in_blocks:
					snps_in_blocks.append(j)
	snps_in_blocks.sort()
	#If the filter option is empty or if it begins with 'n' for non-recombinant, any SNP in the input list
	#that does not fall within a Gubbins block will be output to file. If the option is 'r' for recombinant,
	#any snp that does fall within a Gubbins block will output to file.
	#If the option does not begin with 'r' or 'n', the script will exit and ask the user to specify option
	if filter == None or filter[0].lower() == 'n':
		with open(prefix+'GubbinsFilteredSNPsNotInRecombBlocks.txt', 'w') as output_handle:
			c=0
			for k in snps:
				if k not in snps_in_blocks:
					output_handle.write(str(k)+"\n")
					c+=1
			filename = str(output_handle)
			fname_list = filename.replace(" ",",").split(",")
			print ""
			print "...", c, "Non-recombinant SNP positions written to file", str(fname_list[2]), "..."
	elif filter[0].lower() == 'r':
		with open(prefix+'GubbinsFilteredSNPsInRecombBlocks.txt', 'w') as output_handle:
			c=0
			for k in snps:
				if k in snps_in_blocks:
					output_handle.write(str(k)+"\n")
					c+=1
			filename = str(output_handle)
			fname_list = filename.replace(" ",",").split(",")
			print ""
			print "...", c, "Recombinant SNP positions written to file", str(fname_list[2]), "..."
	else:
		print "\nError! Filter option not recognised. Re-run script without input for -f switch or specify 'recombinant' (or just 'r', upper- or lower-case) on 'nonrecombinant' (or just 'n', upper- or lower-case) (see help on -h).\n"
		sys.exit(0)


block_ranges = get_blocks(args.blocks_file)
snps_in_list = get_snps(args.snp_list)
filter_snps(block_ranges, snps_in_list, args.filter_option, args.output_filename_prefix)

print "\nDone.\n"
