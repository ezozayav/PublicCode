#!/usr/bin/env python

"""
13th Apr 2016
dr.mark.schultzm@gmail.com
github: schultzm
Find the line containing the best taxonomic hit in Kraken output tables.
"""

import argparse
import os

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Find the best matched taxonomic \
                                 hit in the input kraken table(s).")
PARSER.add_argument("-i", "--input", help="Input kraken table file(s)",
                    nargs="+", required=True)
PARSER.add_argument("-n", "--nlines", help="Number of lines to search infiles",
                    type=int, default=10, required=False)
PARSER.add_argument("-t", "--tax_level", help="D, P, C, O, F, G, S",
                    default='S', required=False)
ARGS = PARSER.parse_args()

def file_head(file, nlines):
    """
    Will grab the first n lines from input file and store it in 'head'.
    """
    with open(file, 'r') as input_handle:
        head = [next(input_handle).rstrip().split('\t') \
                for x in xrange(nlines)]
    return head

def get_taxon(file_head, tax_level):
    """
    Will locate the best hit in the kraken table file 'head'.
    """
    #Taxonomic classifications inverted.
    invert_hierarchy = ["S", "G", "F", "O", "C", "P", "D"]
    hit = []
    for i in invert_hierarchy:
        if len(hit) == 0:
            for j in file_head:
                if i in j:
                    hit.append(j)
    #strip leading whitespaces from elements in list.
    hit = [[s.lstrip() for s in inner] for inner in hit]
    return hit

for i in ARGS.input:
    best_hit = get_taxon(file_head(i, ARGS.nlines), ARGS.tax_level)
    #Flatten list by one level and print it
    best_hit = '\t'.join([str(x) for x in best_hit])
    #Print result to screen.  If desired, output to file with '> output.txt'.
    print i, best_hit
