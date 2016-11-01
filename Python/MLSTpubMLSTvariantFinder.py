#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
dr.mark.schultz@gmail.com
31-Oct-2016
github: schultzm

Given a table of ST allelic profiles (e.g., from pubMLST),
output any ST with <= this many loci differing from the reference.
Example ST table here:
wget https://goo.gl/H8Y6CZ -O EcloSTprofiles.tsv


Example usage:
    python MLSTpubMLSTvariantFinder.py -s 114 -d 1 -t EcloSTprofiles.tsv

'''

import argparse
import pandas as pd

#set up the argument parser
PARSER = argparse.ArgumentParser(description="Finds all variants within a \
                                 designated range of a reference ST")
PARSER.add_argument('-s', '--st_number', type=int, help='ST', required=True)
PARSER.add_argument('-v', '--variant_distance', type=int, help='Up to this \
                    many of loci will be different to the reference alleles.',
                    required=True)
PARSER.add_argument('-t', '--table_of_STs', help='Space-, comma- or \
                    tab-delimited text file from pubMLST.  \
                    Expects a header row.', required=True)
ARGS = PARSER.parse_args()

def zip_rows(dframe, refst, nonrefst, variant_dist):
    '''
    Zip two rows of a pandas dataframe.
    For each tuple, return 1 if refst != nonrefst. Sum the
    list.  Return True if sum is less than or equal
    to variant_dist.
    '''
    if sum([1 if i[0] != i[1] else 0 for i in zip(dframe.loc[refst],
                                                  dframe.loc[nonrefst])][:-1]
                                                  ) <= variant_dist:
        return True


def find_variants(st_file):
    '''
    Compare the allelic profile of each ST to the reference ST, including
    the reference ST to itself.
    '''
    #store file as a pandas df
    st_tab = pd.read_csv(st_file, header=0, index_col=0, sep=None,
                         engine='python')
    #rows to keep
    rows_within_variant_dist = []
    for i in range(0, len(st_tab.index)+1):
        if i in st_tab.index:
            #run the zip_rows function
            if zip_rows(st_tab, ARGS.st_number, i, ARGS.variant_distance):
                rows_within_variant_dist.append(i)
    #Keep the rows with True
    return st_tab.loc[rows_within_variant_dist, :]


def main():
    '''
    Read in the file. Find the rows within the range. Print the rows as tsv.
    '''
    table = find_variants(ARGS.table_of_STs)
    print table.to_csv(sep='\t')


if __name__ == '__main__':
    main()
