#!/usr/bin/env python

'''
Concatenate multiple nexus alignments into a 'charset'-defined supermatrix.
dr.mark.schultz@gmail.com
github schultzm
20160624_1220
'''

import argparse
import os
import sys
from Bio import AlignIO
from Bio.Nexus import Nexus
# from Bio.Alphabet import generic_dna

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Concatenate multiple nexus \
                                 alignments.  Handles many input alignments \
                                 (with e.g., *.nexus).')
PARSER.add_argument('-a', '--alignments', help='Input alignments.', nargs='+',
                    required=True)
PARSER.add_argument('-o', '--outname', help='Output filename', default=None, \
                    required=False)
ARGS = PARSER.parse_args()

def concatenate(file_list, outname):
    nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
    combined = Nexus.combine(nexi)
    if outname == None:
        combined.write_nexus_data(filename=open('supermatrix.nex', 'w'))
    else:
        combined.write_nexus_data(filename=open(outname, 'w'))

concatenate(ARGS.alignments, ARGS.outname)
