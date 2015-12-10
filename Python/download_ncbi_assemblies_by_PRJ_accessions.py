##########################
#dr.mark.schultz@gmail.com
#github schultzm
#date 06/11/15
#requires Entrez Direct from http://www.ncbi.nlm.nih.gov/books/NBK179288/
##########################

import os
import argparse

#set up the argument parser
parser = argparse.ArgumentParser(description = "Will download assemblies from NCBI WGS given accession number(s)")
parser.add_argument('-a', '--accession_ids', nargs ='+', help = 'Accession numbers.  If more than one accession, separate list members by whitespace', required=True)
args = parser.parse_args()

#This function does all the work
def get_assembly(accession):
	"""
	Explaining the pipe stored under the 'command' variable below:
	-using Entrez Direct, search for the accession on NCBI (esearch)
	-grab the sequin record (efetch)
	-parse the sequin record in stdout with grep, searching for the PRJ
	accession
	-remove all the junk around the grepped text (using sed to get rid of the
	braces, spaces and inverted commas)
	-pipe this PRJ number to a second round of searching NCBI (using xargs and
	esearch)
	-fetch the fasta records associated with the PRJ accession and store them
	in an mfasta file

	For an explanation on the way accessions and PRJ numbers are linked to assemblies, see:
	http://www.ncbi.nlm.nih.gov/genbank/wgs
	"""
	command = "esearch -db nucleotide -query "+accession+" | efetch --format runinfo | grep 'PRJ' | sed -E 's/\"//g;s/}//g;s/,//g;s/ //g' | xargs esearch -db nucleotide -query | efetch --format fasta > "+accession+".mfasta"
	print "downloading assembly (mfasta format) for accession "+accession+" ..."
	os.system(command)

#execute the get_assembly function for each accession
for i in args.accession_ids:
	get_assembly(i)

LENGTH = len(args.accession_ids)
print '\nDone. '+str(LENGTH)+' files written.\n'

