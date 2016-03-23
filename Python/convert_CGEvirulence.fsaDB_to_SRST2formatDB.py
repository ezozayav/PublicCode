#!/usr/bin/env python

"""
This script will convert CGE virulence gene databases to SRST2 compatible
format.
Carefully check the input .fsa file from CGE to check line breaks.  A quick
way to check and get the line number is to do
"grep '>' *.fsa -n".  Then scroll through the screen output and check for
strange line breaks.  Fix accordingly, then run this script.

gtihub: schultzm
email: dr.mark.schultz@gmail.com
date: 23rd March 2016
"""

import argparse
from itertools import izip

PARSER = argparse.ArgumentParser(description="Convert CGE database to \
                                 SRST2 db.")
PARSER.add_argument("-i", "--input", help="Input databases (fasta format)",
                    nargs='+', required=True)
ARGS = PARSER.parse_args()

def rearrange_header(infile):
    gene = []
    with open(infile, "r") as input_handle:
        for line in input_handle:
            if '>' in line:
                line = line.rstrip().replace('>','>:').split(':')
                gene.append(line[1])
    gene = list(set(gene))
    gene = sorted(gene, key = str.lower)
    gene_dict = {}
    c=1
    for i in gene:
        gene_dict[i]=str(c)
        c+=1
    with open(infile, "r") as input_handle:
        with open(infile+"_SRST2format.fsa", "w") as output_handle:
            for line in input_handle:
                if '>' in line:
                    line = line.rstrip().replace('>','>:').split(':')
                    line_rearranged=line[0]+gene_dict[line[1]]+'__'+line[1]+\
                                    '__'+line[1]+'_'+line[2]+'__'+\
                                    ':'.join(line[3:])+'\n'
                    output_handle.write(line_rearranged)
                else:
                    line=line.rstrip()
                    output_handle.write(line+'\n')

for i in ARGS.input:
    rearrange_header(i)
