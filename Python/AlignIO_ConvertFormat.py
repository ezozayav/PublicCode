#!/usr/bin/env python

'''
Convert DNA alignments between BioPython-supported formats.
dr.mark.schultz@gmail.com
github schultzm
20160624_1120
'''

import argparse
import os
import sys
from Bio import AlignIO
from Bio.Alphabet import generic_dna

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Convert DNA alignments between \
                                 BioPython-supported formats.  Handles \
                                 many input alignments (with e.g., *.fasta).')
PARSER.add_argument('-a', '--alignment', help='Input alignment(s)', nargs='+',
                    required=True)
PARSER.add_argument('-i', '--informat', help='Alignment input format (default \
                    fasta)', default='fasta', required=False)
PARSER.add_argument('-o', '--outformat', help='Alignment output format \
                    (default fasta)', default='fasta', required=False)

ARGS = PARSER.parse_args()

def read_alignment(alignment, informat, outformat):
    align = AlignIO.read(alignment, informat, alphabet=generic_dna)
    out_basename = os.path.splitext(alignment)[0]
    outname = out_basename+'.'+outformat
    if os.path.exists(outname):
        print outname+' already exists. Nothing to do. Moving to next request.'
    if outname != alignment:
        with open(outname, 'w') as output_handle:
            AlignIO.write(align, output_handle, outformat) 
            print '\nConverted '+informat+' to '+outformat
            print 'Written file to '+outname+'. Here is a preview:'
            print align
            print ''
for i in ARGS.alignment:
    read_alignment(i, ARGS.informat, ARGS.outformat)
