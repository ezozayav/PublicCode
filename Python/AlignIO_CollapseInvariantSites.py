###########################
#dr.mark.schultz@gmail.com#
#120114####################
###########################

#import modules
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus
from Bio.Seq import Seq
import argparse
import glob
import itertools


#parse command line input
parser = argparse.ArgumentParser(description = "Will read in an alignment in any format, and spit out only the variant sites to a new nexus-formatted alignment.")
parser.add_argument('-n', '--filename', help = "Name of input file.", required = True)
parser.add_argument('-i', '--informat', help = "File format of input. e.g., genbank, fasta, phylip", required = True)
parser.add_argument('-g', '--gap_char', help = "Exact gap character being used. Default = 'None'", default=None, required = False)
args = parser.parse_args()

#convert files
def read_collapse(file, informat, gapchar):
	with open(file, 'r') as input_handle:

		alignment = AlignIO.read(input_handle, informat, alphabet=generic_dna)
		summary_align = AlignInfo.SummaryInfo(alignment)
		first_seq = (alignment[0].seq)
		length_first_seq = len(first_seq)

# 		chars_to_ignore = ['N']
		my_pssm = summary_align.pos_specific_score_matrix(first_seq)

		index = 0
		count = 0
		invariant_sites_counter = 0
		invariant_position_index = []

		for i in my_pssm.pssm:
			A = i[1]['A']
			C = i[1]['C']
			G = i[1]['G']
			T = i[1]['T']
			if gapchar != None:
				print gapchar
				gap = i[1][gapchar]
				x = [gap, A, C, G, T]
			if gapchar == None:
				x = [A, C, G, T]
				print x
			y = []
			for j in x:
				if j > 0:
					y.append(1)
				else:
					y.append(0)
			if sum(y[1:len(y)]) > 1:
				pass
			else:
				invariant_sites_counter += 1
				invariant_position_index.append(count)
			count += 1

		alignment_indices_to_write = []
		n_alignments = []

		for i in range(0,length_first_seq):
			if i not in invariant_position_index:
				alignment_indices_to_write.append(i)

		def ranges(i):
			for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
				b = list(b)
				yield b[0][1], b[-1][1]

		blocks = list(ranges(alignment_indices_to_write))
		print '\nExcluding', str(len(invariant_position_index)),'sites at positions:',invariant_position_index,'\n'
		print 'Including sites at positions:',blocks,'\n'
		for i in blocks:
			alignment_iteration = MultipleSeqAlignment(alignment[:,i[0]:i[1]+1], alphabet = generic_dna).format('nexus')
			n_alignments.append(('site'+str(i[0])+'to'+str(i[1]+1),Nexus.Nexus(alignment_iteration))) #

		#combine the alignments in n_alignments
		combined = Nexus.combine(n_alignments)
		with open(file+'_collapsed.nexus', 'w') as output_handle:
			print 'Writing collapsed alignment to:',file+'_collapsed.nexus\n'
			combined.write_nexus_data(output_handle)

read_collapse(args.filename, args.informat.replace('.',''), args.gap_char)

print "\nDone.\n"
