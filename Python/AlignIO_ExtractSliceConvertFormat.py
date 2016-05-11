#!/usr/bin/env python

"""
Extract sub-alignment from master alignment by alignment position (column) 
range (start=x1, stop=x2).  Converts between alignment formats.
"""

import argparse
from Bio import AlignIO
from Bio.Alphabet import generic_dna
import os
import sys

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Extract sub-alignment from master \
                                 alignment by alignment position (column) \
                                 range (start=x1, stop=x2). Converts between \
                                 formats.")
PARSER.add_argument("-a", "--alignment", help="Input alignment (BioPython \
                    supported formats)",
                    required=True)
PARSER.add_argument("-i", "--informat", help="Alignment input format (default \
                    fasta)", default="fasta", required=False)
PARSER.add_argument("-o", "--outformat", help="Alignment output format \
                    (default fasta)", default="fasta", required=False)
PARSER.add_argument("-b", "--begin_position", help="Start slice", 
                    type=int, required=True)
PARSER.add_argument("-e", "--end_position", help="End slice", 
                    type=int, required=True)

ARGS = PARSER.parse_args()

def read_alignment(alignment, informat, outformat, start, stop):
    align = AlignIO.read(alignment, informat, alphabet=generic_dna)
    out_basename = os.path.splitext(alignment)[0]
    algn_length = align.get_alignment_length()
    print "\nInput alignment is "+str(algn_length)+" characters."
    end_pos = stop
    if stop>algn_length:
        print "\nNB: you have requested an end position beyond the "+\
               "length of the alignment.  "
        end_pos = algn_length
    if stop<start or start<0:
        print "\nFatal: your begin and end positions need re-assessment."+\
              "  Exiting now."
        print ""
        sys.exit()
    outname = out_basename+"_pos"+str(start)+"to"+str(end_pos)+"."+outformat
    with open(outname, "w") as output_handle:
        algn = align[:, start:stop]
        AlignIO.write(algn, output_handle, outformat) 
        print "\nExtracted "+outformat+"-formatted sub-alignment from "+\
        "positions "+str(start)+" to "+str(end_pos)+" and written it to "+\
        outname+".  Here is a preview:"
        print ""
        print algn
        print ""

read_alignment(ARGS.alignment, ARGS.informat, ARGS.outformat, 
               ARGS.begin_position, ARGS.end_position)
