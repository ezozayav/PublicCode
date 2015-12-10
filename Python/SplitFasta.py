########
#dr.mark.schultz@gmail.com
#080114
########

#import the modules
from Bio import SeqIO
import argparse
import glob

#set up the arguments parser to deal with the command line input and provide help text
parser = argparse.ArgumentParser(description = "Crude script to read in all files with \'.fasta\' (must be exactly that) extension and split into individual output files containing each fasta in each input file.  Will only read .fasta files from the `pwd`.")
args = parser.parse_args()

filelist = glob.glob('*.fasta')

for i in filelist:
	with open(i, 'r') as input_handle:
		for seq_record in SeqIO.parse(input_handle, 'fasta'):
			with open(i.replace('.fasta','')+'_'+seq_record.id+'.fasta', 'w') as output_handle:
				SeqIO.write(seq_record, output_handle, 'fasta')


