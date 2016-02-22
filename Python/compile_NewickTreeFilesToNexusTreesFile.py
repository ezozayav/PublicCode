#!/usr/bin/env python

"""
22nd Feb 2016
dr.mark.schultzm@gmail.com
github: schultzm

Input tree files containing single trees to receive a compiled file \
containing all input trees.
"""

import argparse
import os
import csv
import linecache

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Compile multiple newick tree \
                                 files into a single nexus file")
PARSER.add_argument("-i", "--input", help="Input trees",
                    nargs="+", required=True)
ARGS = PARSER.parse_args()

def compile(tree_files):
    with open("TreesCompiled.nex.trees", "w") as output_handle:
        output_handle.write("#NEXUS\nbegin trees;\n")
        for tree_file in tree_files:
            tree_file_prefix = os.path.splitext(os.path.basename(tree_file))[0]
            tree_name = tree_file_prefix.replace(".","_").split("_")[8]
            tree = linecache.getline(tree_file, 1)
            output_handle.write("tree "+tree_name+" = [&U] "+str(tree))
        output_handle.write("End;")

compile(ARGS.input)
