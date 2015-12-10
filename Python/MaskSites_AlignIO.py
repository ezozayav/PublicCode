###########################
#dr.mark.schultz@gmail.com#
#090415####################
###########################

#import modules
from Bio import AlignIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import argparse
import os

#parse command line input
parser = argparse.ArgumentParser(description = "Will read in an alignment, will read a blocks file, will mask positions in the blocks file with 'NNNN'.")
parser.add_argument('-a', '--alignment_file', help = "Name of input alignment file.", required = True)
parser.add_argument('-i', '--informat', help = "File format of input. e.g., genbank, fasta, phylip", required = True)
parser.add_argument('-o', '--outformat', help = "File format of output. e.g., fasta, genbank, phylip", required = True)
parser.add_argument('-b', '--blocks_file', help = "Block coordinates, named 'xxxx.yyyy.importation_status.txt', with first column strain name, second column block beginning, third column block end")
args = parser.parse_args()

#here we will set a global variable containing keys as taxon names, and values as a list of block position tuples
dict_of_blocks = {}

#will write the block positions to dict_of_blocks
def store_blocks(blocks):
	with open(blocks, 'r') as blocks_handle:
		next(blocks_handle)
		for line in blocks_handle:
			row_vals = [i for i in line.rstrip('\n').split('\t')]
			if row_vals[0] not in dict_of_blocks:
				dict_of_blocks[row_vals[0]]=list()
			dict_of_blocks[row_vals[0]].append((int(row_vals[1]), int(row_vals[2])))

#convert files
def alignment_mask(alignment_file, informat, outformat, blocks):
	#next line must be AlignIO.read and not AlignIO.parse
	alignment = AlignIO.read(open(alignment_file), informat)
	#iterate through the record.seq objects and convert
	print ""
	for record in alignment:
		#convert the sequence from immutable to mutable
		record.seq = record.seq.tomutable()
		for key, value in dict_of_blocks.iteritems():
			if key == record.id:
				blocks_to_mask = value
				print "In sequence", record.id, "masking blocks:", blocks_to_mask, "\n"
				for m in blocks_to_mask:
					record.seq[m[0]:m[1]] = (m[1]-m[0])*'?'
	with open(os.path.splitext(alignment_file)[0]+"_masked."+outformat, 'w') as output_handle:
		print "Writing masked alignment (recombinant values masked with 'N') to", output_handle, "\n"
		AlignIO.write(alignment, output_handle, outformat)

alignment_mask(args.alignment_file, args.informat.replace('.',''), args.outformat.replace('.',''), store_blocks(args.blocks_file))

print "\nDone.\n"
