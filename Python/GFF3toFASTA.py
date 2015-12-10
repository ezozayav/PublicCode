import argparse

#set up the arguments
parser = argparse.ArgumentParser(description = "Converts GFF3 files to multi-fasta")
parser.add_argument('-a', '--accession_ids', nargs ='+', help = 'Accession numbers.  If more than one accession, separate list members by whitespace', required=True)
args = parser.parse_args()

from BCBio import GFF

def convert_gff(accession):
	in_file = accession

	in_handle = open(in_file)
	with open(in_file, 'r') as input_handle:
		infile_name = in_file.replace('.gff','')
		print "\nProcessing "+in_file+"..."
		with open(infile_name+'.mfasta', 'w') as output_handle:
			print '\tExtracting contigs to '+infile_name+'.mfasta'
			for rec in GFF.parse(input_handle):
				output_handle.write('>'+str(infile_name)+'_'+str(rec.id)+'\n'+str(rec.seq)+'\n')
				print '\t\tWritten '+str(rec.id)

for i in args.accession_ids:
	convert_gff(i)

print 'Done\n'
