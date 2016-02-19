#!/usr/bin/env python

"""
19th Feb 2016
dr.mark.schultzm@gmail.com
github: schultzm

Input one or more alleles.csv files (with sorted SNP positions), convert to
fasta files.
"""

import argparse
import os
import csv
# import linecache

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Convert alleles.csv file(s) \
                                 to fasta format")
PARSER.add_argument("-i", "--input", help="Input alleles.csv file(s)",
                    nargs="+", required=True)
ARGS = PARSER.parse_args()

def convert(alleles_csv):
    """
    Convert an alleles.csv file (e.g., from RedDog mapping pipe) to
    .fasta format file
    """
    alleles_data = list(csv.reader(open(alleles_csv, "r")))
    alleles_data_zipped = zip(*alleles_data)
    input_handle_prefix = os.path.splitext(os.path.basename(alleles_csv))[0]
    with open(input_handle_prefix+".fasta", "w") as output_handle:
        for i in alleles_data_zipped[1:]:
            output_handle.write(">"+i[0]+"\n"+''.join(i[1:])+"\n")

for ALLELE_FILE in ARGS.input:
    convert(ALLELE_FILE)
