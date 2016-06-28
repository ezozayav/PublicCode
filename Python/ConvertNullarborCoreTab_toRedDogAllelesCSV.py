#!/usr/bin/env python

"""
20160628
dr.mark.schultzm@gmail.com
github: schultzm
Converts tseemann Nullarbor formatted core.tab file to dedwards RedDog 
Alleles.csv format.  Outputs one file per alignment contig.
"""

import argparse
import os
import random

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Convert a Nullarbor formatted \
                                 core.tab file to RedDog alleles.csv \
                                 formatted file.")
PARSER.add_argument("-i", "--input", help="Input core.tab file (s)",
                    nargs="+", required=True)
PARSER.add_argument("-r", "--remove_reference", help="Would you like to remove \
                    the reference sequence from the table [Y, N] ?  NB: this \
                    may leave you with invariant sites. Default 'N'.  Assumes \
                    reference is in the third (1-based index) column.", \
                    default='n', required=False)
ARGS = PARSER.parse_args()

def get_contig_names(core_tab_file):
    """
    Determines the names of each contig in the core_tab_file and returns them.
    """
    with open(core_tab_file, "r") as input_handle:
        next(input_handle)
        contigs = [line.rstrip("\n").split("\t")[0] for line in input_handle]
        file_sets = set(contigs)
    return file_sets

def split_core_tab(core_tab_file, file_contigs, remove_ref):
    """
    Writes out an alleles.csv for each contig in the core_tab_file.
    """
    base = os.path.splitext(core_tab_file)[0]
    for i in file_contigs:
        with open(core_tab_file, "r") as input_handle:
            first = next(input_handle)
            first = first.rstrip("\n").split("\t")[1:-4]
            if 'y' in remove_ref.lower():
                first.pop(1)
            illegal_chars = '\/:*?"<>|^ '
            cntig_name = i
            for char in illegal_chars:
                if char in cntig_name:
                    cntig_name = cntig_name.replace(char, "_")
            with open(base+"_"+cntig_name+"_alleles.csv", "w") as output_handle:
                output_handle.write(",".join(first)+"\n")
                for line in input_handle:
                    line_vals = line.rstrip("\n").split("\t")
                    if line_vals[0] == i:
                        if 'y' in remove_ref.lower():
                            line_vals.pop(2)
                        output_handle.write(",".join(line_vals[1:-4])+"\n")
    if 'y' in remove_ref.lower():
        print '\nRemove reference requested.  You may need to remove ' \
              'invariant sites.'

for file in ARGS.input:
    split_core_tab(file, get_contig_names(file), ARGS.remove_reference)

#End message
int_val = random.randint(1,5)
mssg1 = "\nDone. Thanks for using.\n"
mssg2 = "\nAll done.\n"
mssg3 = "\nDone. Check your working directory for the output files.\n"
mssg4 = "\nDone. Email dr.mark.schultz@gmail.com for complaints.\n"
mssg5 = "\nDone. Now, what does it all mean?\n"

if int_val == 1:
    print mssg1
if int_val == 2:
    print mssg2
if int_val == 3:
    print mssg3
if int_val == 4:
    print mssg4
if int_val == 5:
    print mssg5
