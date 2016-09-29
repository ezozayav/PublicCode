#!/usr/bin/env python


'''
Query entrez direct for full scientific rank containing the query.
Dependencies: install Entrez direct
(https://www.ncbi.nlm.nih.gov/books/NBK179288/).

Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20160915_1300
'''


import argparse
from multiprocessing import cpu_count, Pool
from subprocess import Popen, PIPE
import shlex


VERSION = '0.0.1'


# set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Get full NCBI scientific name.')
PARSER.add_argument('-i', '--infile', help="File of taxon names, \
                    one per line", required=True)
ARGS = PARSER.parse_args()

def query_entrez(taxon):
    '''
    Given a taxon name, query NCBI entrez direct for the full scientific rank
    for that taxon.
    '''
    cmd = "esearch -db taxonomy -query '"+taxon+"'"
    cmd2 = "efetch -format xml"
    cmd3 = "xtract -pattern Taxon -element ScientificName"
    args = shlex.split(cmd)
    args2 = shlex.split(cmd2)
    args3 = shlex.split(cmd3)
    proc1 = Popen(args, stdout=PIPE)
    proc2 = Popen(args2, stdin=proc1.stdout, stdout=PIPE)
    proc3 = Popen(args3, stdin=proc2.stdout, stdout=PIPE)
    output = proc3.stdout.read()
    return output.rstrip().split('\t')

def read_infile(infile):
    '''
    Read the lines of the infile. Store each line.
    '''
    taxa = filter(None, [line.rstrip() for line in open(infile).readlines()])
    return taxa


def main():
    '''
    In parallel, execute the query_entrez function.
    '''
    taxon_list = read_infile(ARGS.infile)
    n_cpus = cpu_count()
    if len(taxon_list) < n_cpus:
        n_procs = len(taxon_list)
    else:
        n_procs = n_cpus
    proc_pool = Pool(n_procs)
    results = proc_pool.map(query_entrez, taxon_list)
    for i in results:
        print '\t'.join(i)


if __name__ == '__main__':
    main()
