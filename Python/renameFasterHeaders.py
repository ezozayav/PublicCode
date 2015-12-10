#First, import the modules
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import glob

#set up the arguments parser to deal with the command line input
parser = argparse.ArgumentParser(description = "replaces fasta header inside the file with that file's filename")
parser.add_argument('-e', '--filesuffix', help = "e.g., '.fasta'", required = True)

args = parser.parse_args()

def readconvert(filesuffix):
	filelist = glob.glob('*.'+str(filesuffix.replace('.','')))
	for i in filelist:
		with open(i, 'r') as input_handle:
			#NB: for genbank formatted output, must add ", generic_dna" below.
			for seq_record in SeqIO.parse(input_handle, 'fasta', generic_dna):
				with open(i+'rename.fasta', 'w') as output_handle:
					output_handle.write(str('>'+i.replace('.fasta','')+'\n'))
					output_handle.write(str(seq_record.seq)+'\n')

readconvert(args.filesuffix)
