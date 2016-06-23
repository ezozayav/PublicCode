#!/usr/bin/env python

'''
Extract CDS DNA sequences from gbk files.
20160623_1120
github schultzm
dr.mark.schultz@gmail.com
'''

import argparse
from Bio.Alphabet import generic_dna
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Get CDS DNA seq(s) from genbanks')
PARSER.add_argument('-c', '--CDS_names', help='Names of CDSs', nargs='+',
                    required=False)
PARSER.add_argument('-g', '--genbank_files', 
                    help="Input genbanks sequence(s).", nargs='+', 
                    required=True)
ARGS = PARSER.parse_args()

def getDNAseq(cds_name, gbk):
    '''
    Reads in the requested CDS name, searches for the CDS in the genbank file, 
    returns the DNA sequence, with isolate name and product name in the 
    sequence description.
    '''
    basename = os.path.splitext(os.path.basename(gbk))[0]
    n_hits = 0
    for gb_record in SeqIO.parse(open(gbk,'r'), 'genbank'):
        for index, feature in enumerate(gb_record.features):
            if feature.type == 'CDS':
                gb_feature = gb_record.features[index]
                if 'gene' in gb_feature.qualifiers:
                    if cds_name == gb_feature.qualifiers['gene'][0]:
                        product = gb_feature.qualifiers['product']
                        DNAseq = gb_feature.extract(gb_record.seq)
                        record_new = SeqRecord(Seq(str(DNAseq), generic_dna), \
                                     id=basename, name=cds_name, \
                                     description=product[0])
                        return record_new
                        n_hits += 1
    #Tell the user if copy number of the CDS is greater than single-copy
    if n_hits > 1:
        print 'Warning, '+n_hits+' found in '+gbk+' for '+cds_name+'.'

def writeFasta(seqs, cds_name):
    '''
    Takes the list of sequences stored in SEQS global and writes them to file.
    '''
    with open(cds_name+'.fasta', 'w') as output_handle:
        for i in seqs:
            SeqIO.write(i, output_handle, 'fasta')

#For every locus, check each genbank file and return the locus.
#Write the locus to file.
for i in ARGS.CDS_names:
    SEQS = []
    for j in ARGS.genbank_files:
        print j
        seq = getDNAseq(i, j)
        SEQS.append(seq)
    #use filter on SEQS to get rid of 'None' objects in list
    writeFasta(SEQS, i)

print '\nDone\n'
