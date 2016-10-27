#!/usr/bin/env python


'''
Extract a list of isolates from a BioPython supported alignment file.
Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20160715_1200
'''


import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna


#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Extract subset of taxa from an \
                                 alignment file (BioPython supported \
                                 formats).  Include or exclude the taxa.')
PARSER.add_argument('-a', '--alignment_file', help='Feed me the alignment!',
                    required=True)
PARSER.add_argument('-i', '--informat', help='BioPython supported \
                    formats (use whole words).', required=True)
PARSER.add_argument('-t', '--taxon_list', help='Single column text file \
                    containing the taxon IDs to subset.', required=True)
PARSER.add_argument('-x', '--include_or_exclude', help='Do you want to \
                    \'include.\' or \'exclude\' the taxa?', required=True)
PARSER.add_argument('-o', '--outformat', help='Name of outfile.',
                    required=True)
ARGS = PARSER.parse_args()


def main():
    '''
    Read in the taxa.  Read in the alignment.  Subset records in the \
    alignment if they match or don't match the taxa.  Export the subset \
    alignment in requested outformat.
    '''
    taxa = read_taxa(ARGS.taxon_list)
    aln = read_alignment(ARGS.alignment_file, ARGS.informat.lower())
    sub = subset(aln, taxa, ARGS.include_or_exclude)
    print sub
    outfile_name = ARGS.alignment_file+'_'+str(len(sub))+'isolate_subset.'+\
               ARGS.outformat.lower()
    with open(outfile_name, 'w') as outfile:
        AlignIO.write(sub, outfile, ARGS.outformat.lower())
        print '\nAlignment written to '+outfile_name+'\n'


def read_taxa(t_list):
    'Return the list of taxa from the taxon_list file.'
    return filter(None, [line.rstrip() for line in open(t_list).readlines()])


def read_alignment(infile, informat):
    'Return the alignment.'
    return AlignIO.read(open(infile), informat, alphabet=generic_dna)


def subset(alignment, taxa, include_or_exclude):
    'Extract the subset of included or excluded taxa'
    if include_or_exclude.lower() == 'include':
        sub_set = [record for record in alignment if record.id in taxa]
        print '\nIncluding isolates in taxon_list infile...\n'
    if include_or_exclude.lower() == 'exclude':
        sub_set = [record for record in alignment if record.id not in taxa]
        print '\nExcluding isolates in taxon_list infile...\n'
    return MultipleSeqAlignment(sub_set)

if __name__ == '__main__':
    main()
