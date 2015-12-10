#Sep 23 2014
#dr.mark.schultz@gmail.com
#github: schultzm
#This file will split a SNP alignment into overlapping user-selected formatted sub-alignments.'


#First, import the modules
import argparse
import os
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus
# print dir(generic_dna)

#sys.exit(0)
parser = argparse.ArgumentParser(description = "This file will split a SNP alignment into overlapping user-selected formatted sub-alignments.  File name will indicate slice position.")
parser.add_argument('-i', '--input', help = 'Input alignment file to pull the sub-alignments from.', required=True)
parser.add_argument('-f', '--informat', help = 'Informat of input file.', required = True)
parser.add_argument('-o', '--outformat', help = 'Outformat of output file (Biopython supported alignments, e.g., fasta, nexus or phylip). NB: if specifying phylip outformat, input-file sequence names will be truncated to <= 10 characters in the output files, which may result in non-unique sequence names; thus, being a problem for downstream analyses.', required = True)
parser.add_argument('-w', '--SNPs_in_window', type = int, help = 'Size of window output to each file (i.e., output files will contain this number of SNPs, up to the remainder alignment, which will contain less).', required = True)
parser.add_argument('-s', '--slide', type = int, help = 'Slide this number of sites to direction 3\' on subsequent sub-alignment extractions.', required = True)

args = parser.parse_args()

#list the files to be converted
os.system('pwd')
print 'Present working directory (see immediately above):'
print 'File to be subsampled (-i): %s' % args.input
print 'Input file format (-f): %s' % args.informat.replace('.','')
print 'Output file format (-o): %s' % args.outformat.replace('.','')
print 'SNPs requested per output file (-w): %s' % args.SNPs_in_window
print 'Slide size (number of SNPs): %s' % args.slide

def alignment_slicer(input, informat, outformat, SNPs, slide):
	alignment =  AlignIO.read(input, informat, alphabet = generic_dna)
	alignment_seq_count = len(alignment)
	first_seq = (alignment[0].seq)
	length_alignment = len(first_seq)
	chars_to_ignore = ['N']
	
	start = 0
	end = start + args.SNPs_in_window
	while end <= length_alignment:
		with open(input+'_site'+str(start)+'to'+str(end)+'.'+outformat, 'w') as output_handle:
			
	# 		print 'start:', start
	# 		print 'end:', end
			alignment_iteration = MultipleSeqAlignment(alignment[:, start:end], alphabet=generic_dna)
			if outformat.lower() == 'nexus':
				n_alignments = []
				alignment_iteration = alignment_iteration.format('nexus')
				n_alignments.append(('site'+str(start)+'to'+str(end),Nexus.Nexus(alignment_iteration)))
				combined = Nexus.combine(n_alignments)
				combined.write_nexus_data(output_handle)
			else:
				AlignIO.write(alignment_iteration, output_handle, outformat)
	# 		print alignment_iteration
			start += args.slide
			end += args.slide
	else:
		with open(input+'_site'+str(start)+'to'+str(length_alignment)+'.'+outformat, 'w') as output_handle:
			n_alignments = []
	# 		print 'now in else loop\n'
	# 		print 'start:', start
	# 		print 'end:', length_alignment
			alignment_iteration = MultipleSeqAlignment(alignment[:, start:length_alignment], alphabet=generic_dna)
			if outformat.lower() == 'nexus':
				n_alignments = []
				alignment_iteration = alignment_iteration.format('nexus')
				n_alignments.append(('site'+str(start)+'to'+str(end),Nexus.Nexus(alignment_iteration)))
				combined = Nexus.combine(n_alignments)
				combined.write_nexus_data(output_handle)
			else:
				AlignIO.write(alignment_iteration, output_handle, outformat)
	# 		print alignment_iteration
		print "\ndone\n"
				
alignment_slicer(args.input, args.informat, args.outformat, args.SNPs_in_window, args.slide)
