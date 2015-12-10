#Sep 22 2014
#dr.mark.schultz@gmail.com
#github: schultzm
#This file will return an alignment pssm from an input alignment in any format.  

#First, import the modules
import argparse
import os
import sys
import itertools

from Bio import AlignIO
from Bio.Alphabet import generic_dna
# from Bio import Alphabet
from Bio.Align import AlignInfo
# from Bio.Alphabet import Gapped, Alphabet, generic_dna #IUPAC, 
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus
from Bio.Seq import Seq
# alphabet = generic_dna#Gapped(IUPAC.ambiguous_dna)

# sys.exit(0)
parser = argparse.ArgumentParser(description = "This file will return an alignment Position Specific Score Matrix (PSSM) from an input alignment in any format, will not count gaps as a fifth state.")
parser.add_argument('-i', '--input', help = 'Input alignment file to scan.', required=True)
parser.add_argument('-f', '--informat', help = 'Informat of input file.', required = True)

args = parser.parse_args()

def alignment_input(input, format):
	alignment_pssm =  AlignIO.read(input, format, alphabet = generic_dna)
	alignment_seq_count = len(alignment_pssm)
	first_seq = alignment_pssm[0].seq
	length_first_seq = len(first_seq)
	chars_to_ignore = ['N']
	summary_align = AlignInfo.SummaryInfo(alignment_pssm)
	my_pssm = summary_align.pos_specific_score_matrix(first_seq, chars_to_ignore)
	with open(str(input)+'_PSSM.txt', 'w') as output_handle:
		output_handle.write(str(my_pssm))
		
alignment_input(args.input, args.informat)
