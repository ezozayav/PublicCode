###########################
#dr.mark.schultz@gmail.com#
#080114####################
###########################

#import modules
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import textwrap
import argparse
import glob

#parse command line input
parser = argparse.ArgumentParser(description = "This script will read in all files with specified extension. For each file, it will then spit out the individual sequences as individual files.  Input and output formats can be specified to, for example, convert a multi-fasta file to individual genbank files. File names of outputs will reflect the name of the input files.")
parser.add_argument('-e', '--file_extension', help = "Extension of input file.", required = True)
parser.add_argument('-i', '--informat', help = "File format of input. e.g., genbank, fasta, phylip", required = True)
parser.add_argument('-o', '--outformat', help = "File format of output. e.g., fasta, genbank, phylip", required = True)
args = parser.parse_args()

#get the list of files
def get_file_list(file_extension):
	filelist = glob.glob('*.'+args.file_extension.replace('.',''))
	if len(filelist) == 0:
		print "No input files found..."
		print "...Exiting."
	else:
		print "\nParsing", str(len(filelist)), "input", file_extension.replace('.',''), "files and converting to", args.outformat.replace('.',''), "format."
	return filelist

#convert files
def readconvert(files, informat, outformat):
	for i in files:
		with open(i, 'r') as input_handle:
			#NB: for genbank formatted output, must add ", generic_dna" below.
			for seq_record in SeqIO.parse(input_handle, informat, generic_dna):
				with open(i+'_'+seq_record.id+'.'+outformat, 'w') as output_handle:
					SeqIO.write(seq_record, output_handle, outformat)

#get the file list by calling the 'get_file_list' function and pass the returned output to the readconvert function
files = get_file_list(args.file_extension)
readconvert(files, args.informat.replace('.',''), args.outformat.replace('.',''))

print "\nDone.\n"
