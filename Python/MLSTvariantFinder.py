# -*- coding: utf-8 -*-
"""
#############################
#dr.mark.schultz@gmail.com  #
#08-Jan-2016                #
#github: schultzm           #
#############################

This module takes a reference ST (a combination of alleles) and a table of
isolate STs (allele combinations).  At a user specified distance – that is, the 
number of loci different to the reference ST – the script will output any isolate 
STs that are within this range of the reference.

Example usage:
    python MLSTvariantFinder.py -s 2 2 2 2 2 2 2 -d 4 -t isolate_MLSTs.csv

The above command will execute this module, and compare any isolate MLST
profile in the isolate_MLSTs.csv file to the reference ST of '2 2 2 2 2 2 2'.
The '-d 4' in the command means output to file any isolate within four
differences to the reference ST.  To gather all of the single locus variants, 
use -d 1, the get the doubles use -d 2, etc.
"""

import argparse
import os

#set up the argument parser
PARSER = argparse.ArgumentParser(description="Finds all variants within a \
                                              designated range of a \
                                              reference ST")
PARSER.add_argument('-s', '--reference_sequence_type', nargs='+',
                    help='Reference sequence type (for Pasteur scheme, this \
                          is seven loci, space delimited), \
                          e.g., \'2 2 2 2 2 2 2\'', required=True)
PARSER.add_argument('-d', '--locus_variant_distance',
                    help='Distance (input as an integer = up to this number \
                          of loci different to reference ST) from reference \
                          ST (above) within which to gather other STs from \
                          the ST table (below).', required=True)
PARSER.add_argument('-t', '--table_of_STs',
                    help='Comma-delimited text file with first column as \
                          sample name, second column as ST, then the \
                          subsequent seven columns as the allele number \
                          for each locus. Expects a header row \
                          in input file.', required=True)
ARGS = PARSER.parse_args()

def find_st(st_file):
    """Compare the allele profile of each isolate to the reference ST."""
    with open(st_file, 'r') as input_handle:
        with open(os.path.splitext(st_file)[0]+'_'+ARGS.locus_variant_distance+
                  'LociVariants.tab', 'w') as sample_list:
            # Write the header line of the input file to output file.
            sample_list.write(input_handle.readline())
            # Process all trailing lines in the input file.
            for line in input_handle:
                profile = line.rstrip().split(',')[2:]
                # Compare ST of each sample in input_handle to reference ST
                zipped = zip(profile, ST)
                count = 0
                for i in zipped:
                    if i[0] != i[1]:
                        count += 1
                if count <= int(ARGS.locus_variant_distance):
                    sample_list.write(line)

# Store the reference ST
ST = ARGS.reference_sequence_type

# Execute the find_STs function, comparing the isolate MLSTs to the reference
find_st(ARGS.table_of_STs)

print 'Done'
