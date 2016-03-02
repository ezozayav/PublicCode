#!/usr/bin/env python

"""
02 March 2016
Happy 1st Birthday Mavis! 
dr.mark.schultzm@gmail.com
github: schultzm

Input the alleles csv file (with sorted SNP positions), which has isolate
names in the first row (header).  Subsequent rows contain the sequence data
with each row being a site in the alignment.  The first column of each row is
the SNP position against the reference sequence.  Choose how many SNPs you
would like per tile (i.e., how many rows per file), the script will then copy
the header and that slice to a new file.  Will then move to the next SNP
position after the end of the tile and output that too.  The 'remainder'
alignment may be less than the number of requested SNPs per window.  Also
outputs a file containing a translation between SNP number (as used in the
output file names) and SNP position against the reference.
"""

import argparse
import os
import linecache

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Sample sliding windows along\
                                 the length of an alignment in alleles.csv \
                                 format")
PARSER.add_argument("-i", "--input", help="Input alleles.csv from which to \
                    pull the tiled windows", required=True)
PARSER.add_argument("-w", "--window_size", help="Size of tiled windows \
                    (n SNPs)", type=int, required=True)
ARGS = PARSER.parse_args()

def get_window(alleles_csv, window):
    """
    Main function, gets the tiled windows from the alleles.csv file.
    """
    input_handle_name = os.path.basename(alleles_csv)
    input_handle_prefix = os.path.splitext(input_handle_name)[0]
    window_start = 2
    num_lines = sum(1 for line in open(alleles_csv))-1
    with open("SNPpositionsTranslated_"+input_handle_prefix+".txt", "w") \
    as snp_pos_translated:
        snp_pos_translated.write("SNP_No,SNP_RefPos\n")
        while window_start+window <= num_lines:
            with open(input_handle_prefix+"_SNPnumber"+str(window_start-1)+\
            "to"+str(window_start+window-2)+".csv", "w") as output_handle:
                output_handle.write(linecache.getline(alleles_csv, 1))
                for i in xrange(window_start, window_start+window):
                    line = linecache.getline(alleles_csv, i)
                    snp_pos_translated.write(str(i-1)+\
                    ","+line.rstrip().split(",")[0]+"\n")
                    output_handle.write(line)
            window_start += window
        else:
            with open(input_handle_prefix+"_SNPnumber"+str(window_start-1)+\
            "to"+str(num_lines-1)+".csv", "w") as output_handle:
                output_handle.write(linecache.getline(alleles_csv, 1))
                for i in xrange(window_start, num_lines+1):
                    line = linecache.getline(alleles_csv, i)
                    snp_pos_translated.write(str(i-1)+","+\
                    line.rstrip().split(",")[0]+"\n")
                    output_handle.write(line)

get_window(ARGS.input, ARGS.window_size)
