#!/usr/bin/env python

"""
14th Apr 2016
dr.mark.schultzm@gmail.com
github: schultzm
Converts tseemann Nullarbor formatted core.tab file to dedwards RedDog 
Alleles.csv format.  Outputs one file per 'chromosome'.
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

def split_core_tab(core_tab_file, file_contigs):
    """
    Writes out an alleles.csv for each contig in the core_tab_file.
    """
    base = os.path.splitext(core_tab_file)[0]
    for i in file_contigs:
        with open(core_tab_file, "r") as input_handle:
            first = next(input_handle)
            first = first.rstrip().split("\t")[1:-4]
            illegal_chars = '\/:*?"<>|^ '
            cntig_name = i
            for char in illegal_chars:
                if char in cntig_name:
                    cntig_name = cntig_name.replace(char, "_")
            with open(base+"_"+cntig_name+"_alleles.csv", "w") as output_handle:
                output_handle.write(",".join(first)+"\n")
                for line in input_handle:
                    line_vals = line.rstrip().split("\t")
                    if line_vals[0] == i:
                        output_handle.write(",".join(line_vals[1:])+"\n")

for file in ARGS.input:
    split_core_tab(file, get_contig_names(file))

#End message
int_val = random.randint(1,5)
mssg1 = "Done. Thanks for using."
mssg2 = "All done."
mssg3 = "Done. Check your working directory for the output files."
mssg4 = "Done. Email dr.mark.schultz@gmail.com for complaints."
mssg5 = "Done. Now, what does it all mean?"

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
