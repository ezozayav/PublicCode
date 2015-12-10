##########################
#dr.mark.schultz@gmail.com
#github: schultzm
#170215
##########################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import os
import argparse
import re
import sys, traceback

#example execute "python parse.py -g AAAK01000287.gbk AE016832.gbk > name_of_log_file.txt"
#set up the arguments parser to deal with the command line input
parser = argparse.ArgumentParser(description = "Python/Biopython script for outputting the DNA sequences of all 'genes' in a Genbank input to separate files. Input format Genbank only.  Output format, any biopython supported format (e.g., Genbank, fasta, phylip, nexus).")
parser.add_argument('-g', '--genbank', nargs='+',required=True, help = 'Input names of Genbank file(s), separated by whitespace, to pull the productname and sequence slice from.')
parser.add_argument('-o', '--out_format', default='fasta', required=False, help = "Output file format. Default is 'fasta'.")
args = parser.parse_args()

#parse the genbank file and write every 'gene' in the file out to a separate file with filename and header = Accesion number + gene name + gene location
def extract_genes(file, outformat):
	with open(file, 'r') as handle:
		print "\n********"
		print "Parsing",file,"..."
		accession_no = str(os.path.splitext(file)[0])
		try:
			for record in SeqIO.parse(handle, "genbank"):
			#check if sequence included in genbank file
				if record.seq.count("N")==len(record.seq):
					print "No sequence in file.  Moving onto next file."
					files_without_sequence_data.append(file)
				else:
					alphabet_type = str(record.seq.alphabet)
					if "DNA" in alphabet_type:
						#c will be used to append a number to ensure unique file headers
						c = 1
						#iterate through the record features (e.g., 'gene', 'CDS' etc.)
						for i in range(0, len(record.features)):
							#only consider features that are CDS (coding sequences)
							if record.features[i].type == "CDS":
								#only take information in the CDS annotations from the following qualifiers, if they exist
								#gene name is not always stored in 'gene' annotation, so create a gene name from info stored in any of these keys:
								type_qualifiers = ["gene", "note", "product", "locus_tag", "function"]
								gene_name = []
								for j in type_qualifiers:
									if j in record.features[i].qualifiers:
										#sometimes the feature qualifiers are lists, so convert to single value using 'join'
										#append to the gene_name "container" (list)
										gene_name.append("_".join(record.features[i].qualifiers[j]))
								#setup the regular expression to replace any non-alphanumeric character
								#see https://docs.python.org/3.3/howto/regex.html
								p = re.compile('\W')
								#use re.sub() to replace the non-alphanumeric characters, after joining the gene name list into a single value
								gene_name = p.sub("_", "_".join(gene_name))
								unique_id_prefix = accession_no+"_CDS"+str(c)
								#ensure fasta header is not longer than 30 characters
								fasta_header = unique_id_prefix+"_"+gene_name[0:29-len(unique_id_prefix)]
								#the values stored in the feature variable are needed to get the seq slices
								feature = record.features[i]
								my_seq_location = record.features[i].location
								#for examples on extracting seqs from genbank see http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/#features
								#thanks Peter for the blog
								my_seq = feature.extract(record.seq)
								my_seq_string = str(feature.extract(record.seq))
								record_new = SeqRecord(Seq(my_seq_string, generic_dna), id=fasta_header, description="| 'Original seq coordinates in refseq "+accession_no+": "+str(my_seq_location)+".  Gene ["+gene_name+"]'.")
								outfile_name = fasta_header+"."+outformat
								with open(outfile_name, 'w') as output_handle:
									print "CDS '"+gene_name+"' from accesion "+accession_no+" written to '"+outfile_name+"'."
									SeqIO.write(record_new, outfile_name, outformat)
								c+=1
						successfully_parsed_files.append(file)
					else:
						print "Sequence 'alphabet' is "+alphabet_type+". Unable to extract gene DNA sequence.  Moving onto next file in list."
						files_not_in_DNA_format.append(file)
		#raise error if Genbank file cannot be parsed, and continue without this file.
		except:
			print "ERROR: Genbank parsing of file",file,"failed.  Refer to traceback() error below.  If no error reported below (i.e., you are reading this message as redirected screen output in log file, refer to the traceback() that was output to the screen before opening this file, or re-run this script without redirecting the screen output."
			traceback.print_exc()
			files_failed_parsing.append(file)

#for logging purposes on processing many input files
successfully_parsed_files = []
files_without_sequence_data = []
files_not_in_DNA_format = []
files_failed_parsing = []

print "\nAttempting to process files: "+str(args.genbank)

#execute the extract_genes function on every file in the input list
for i in args.genbank:
	extract_genes(i, args.out_format)

print "\n******** RUN LOG ********\n"

print "Successfully parsed files: ", ", ".join(successfully_parsed_files),"\n"
print "Files without sequence data (no genes to extract): ", ", ".join(files_without_sequence_data),"\n"
print "Files not in DNA format (no DNA sequences to extract): ", ", ".join(files_not_in_DNA_format),"\n"
print "Files failed parsing (e.g., ill-formatted): ", ", ".join(files_failed_parsing),"\n"
print "End.\n"
