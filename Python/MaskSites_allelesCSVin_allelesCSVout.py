###########################
#dr.mark.schultz@gmail.com#
#090915####################
###########################

#import modules
import pandas as pd #deals with the csv. Install from http://pandas.pydata.org/pandas-docs/version/0.15.2/install.html
import argparse
import os

#parse command line input
#example command "python MaskSites_allelesCSVin_FASTAout.py -a x_alleles.csv -b x_alleles_blocks.txt"
parser = argparse.ArgumentParser(description = 'Will read in an "x_alleles.csv" formatted alignment (csv, first column SNP pos as number, trailing columns strains, rowcol values are SNP alleles "A, C, G or T") and a post-ClonalFrameML (see https://code.google.com/p/clonalframeml/) blocks file (named something like "xxxx.yyyy.importation_status.txt"), will mask positions in the alleles.csv with "N" [default] or any other character specified on the "-c" command-line switch.  Strain names must match identically in "blocks_file" and "x_alleles.csv" file.')
parser.add_argument('-a', '--alleles_csv_file', help = "Name of input 'x_alleles.csv' SNP alignment file.", required = True)
parser.add_argument('-b', '--blocks_file', help = 'Block coordinates file, named something like "xxxx.yyyy.importation_status.txt" output from ClonalFrameML, tab-delimited, with first column being strain name, second column block beginning and third column block end.', required=True)
parser.add_argument('-c', '--character_mask', default = 'N', help = 'Character used to mask nucleotides.  If desired, set to e.g., "?", otherwise will default to "N".', required = False)
args = parser.parse_args()

# Set a global variable (type=dictionary) containing keys as taxon names, and
# values as a list containing block-position tuples
dict_of_blocks = {}

# Write the block positions from blocks_file to "dict_of_blocks" global above
def store_blocks(blocks):
	with open(blocks, 'r') as blocks_handle:
		next(blocks_handle)
		for line in blocks_handle:
			row_vals = [i for i in line.rstrip('\n').split('\t')]
			if row_vals[0] not in dict_of_blocks:
				dict_of_blocks[row_vals[0]]=list()
			dict_of_blocks[row_vals[0]].append((int(row_vals[1]), int(row_vals[2])))

#Convert "x_alleles.csv" file to ".fasta" format whilst simultaneously masking SNP positions
def alleles_csv_mask(alleles_csv, blocks, character):
	df = pd.read_csv(alleles_csv)
	colnames = list(df.columns.values)
	#rownames_header=names[0].  This should either be "SNP' or "Pos" in the "x_alleles.csv".
	with open(os.path.splitext(alleles_csv)[0]+'_RecombMasked.fasta', 'w') as output_handle:
		print '\nMasking SNPs in "'+alleles_csv+'" with "'+character+'".'
		print 'Writing masked alignment to "'+output_handle.name+'"...'
		for colname in colnames[1:]:
			#'pos_tuple':  the first value is the SNP position, the second is the allele
			pos_tuple = zip(df[colnames[0]], df[colname])
			if colname in dict_of_blocks:
				blocks_to_mask = dict_of_blocks[colname]
				output_handle.write('>'+colname+'\n')
				for pos in pos_tuple:
					mutate = []
					for block in blocks_to_mask:
						if block[0]<=pos[0]<=block[1]:
							mutate.append('mutate')
					if len(mutate) > 0:
						output_handle.write(character)
					else:
						output_handle.write(pos[1])
				output_handle.write('\n')
			else:
				output_handle.write('>'+colname+'\n')
				for pos in pos_tuple:
					output_handle.write(pos[1])
				output_handle.write('\n')

alignment_mask(args.alignment_file, store_blocks(args.blocks_file), args.character_mask)

print '\nDone.\n'
