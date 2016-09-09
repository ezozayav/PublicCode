#!/usr/bin/env python


'''
Download NCBI files by accession number.
dr.mark.schultz@gmail.com
github schultzm
date 20160909
'''


#import the modules
import argparse
import getpass
import socket
import sys
from Bio import Entrez


#set up the arguments
PARSER = argparse.ArgumentParser(description='Will pull NCBI files' \
                                 ' specified by accession numbers and put' \
                                 ' them into your working directory.')
PARSER.add_argument('-a', '--accession_ids', nargs='+', help='Accession' \
                    ' numbers, white-space delimited', required=True)
PARSER.add_argument('-e', '--user_email', help='User email address',
                    default=None, required=False)
PARSER.add_argument('-f', '--file_format', help='File format. Default' \
                    ' = \'fasta\'.', default='fasta', required=False)
ARGS = PARSER.parse_args()


def get_accession():
    '''
    Given a list of accessions, download the files from NCBI.
    Default format is 'fasta' but will also allow genbank.
    '''
    if ARGS.file_format.lower()[0] == 'g':
        file_type = 'gbwithparts'
        extension = 'gbk'
    else:
        file_type = ARGS.file_format
        extension = ARGS.file_format
    if file_type in {'fasta', 'gbwithparts'}:
        for i in ARGS.accession_ids:
            try:
                handle = Entrez.efetch(db='nucleotide', id=i, rettype=file_type,
                                       retmode='text')
                with open(i+'.'+extension, 'w') as output_file:
                    print 'downloading '+i+' to '+i+'.'+extension
                    output_file.write(handle.read())
            except StopIteration:
                print 'Accession number '+i+' not found'
                continue
            handle.close()
    else:
        print 'File type must be \'genbank\' or \'fasta\'. Exiting now.'
        sys.exit()


def main():
    '''
    Main function.  First tell NCBI who the downloader is via an email address.
    Then run the get_accession function.
    '''
    if ARGS.user_email is None:
        Entrez.email = getpass.getuser()+'@'+socket.getfqdn()
        print '\nDownload owner defaulting to '+Entrez.email
    if ARGS.user_email != None:
        print '\nDownload owner specified as '+ARGS.user_email
        Entrez.email = ARGS.user_email
    get_accession()


if __name__ == '__main__':
    main()
    print 'Done.'
