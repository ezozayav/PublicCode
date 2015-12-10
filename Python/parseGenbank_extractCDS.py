##########################
#dr.mark.schultz@gmail.com
#github: schultzm
#120215
##########################

from Bio import SeqIO
import argparse

#example execute "python parse.py -g AAAK01000287.gbk AE016832.gbk > name_of_log_file.txt"
#set up the arguments parser to deal with the command line input
parser = argparse.ArgumentParser(description = 'This is a script for pulling out CDS to individual files from Genbank')
parser.add_argument('-g', '--genbank', help = 'Input names of Genbank file(s), separated by whitespace, to pull the productname and sequence slice from.', nargs='+',required=True)
args = parser.parse_args()

def split_CDS(file):
	handle = open(file, "rU")
# 	with open('parseGenbank_log.txt', 'w') as log_file:
	for index, record in enumerate(SeqIO.parse(handle, "genbank")):
		c=0
		for i in range(0, len(record.features)):
			if record.features[i].type == 'CDS':
				try:
					#convert the locations in x to an integer list and store in q
					x = str(record.features[i].location).split(":")
					q=[int(j.replace('[','').replace(']','').replace('(+)','').replace('(-)','')) for j in x]
					#get the sequence between the feature.locations in q
					prod = record.features[i].qualifiers['product']
					#fix product names and add unique CDS id to start
					for j in prod:
						#new product name for fasta header = 'k'
						#k will also be the filename
						k = file.replace(".gbk","_")+'CDS'+str(c)+'_'+'pos'+str(q[0])+'to'+str(q[1])+"_"+j.replace('(','_').replace("'","prime").replace(" ","_").replace(")","").replace("-","_").replace(",","_").replace("/","or")
						outfile = str(k)+'.fasta'
						seq = str(record.seq[q[0]-1:q[1]])
						print k+' written to '+outfile
						with open(outfile, 'w') as output_handle:
							#write the seq slice and header to file, manually
							#to do, create a seq object to write 'fasta' format automatically
							output_handle.write('>'+str(k)+'\n'+str(seq)+'\n')
				except:
					#l= error message and file ref
					l='Unable to parse seq locations in CDS'+str(c)+' file '+file
					print l
					continue
			c+=1
	handle.close()

for i in args.genbank:
	split_CDS(i)
