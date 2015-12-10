###########################
#dr.mark.schultz@gmail.com#
#160415####################
###########################
#rename taxa in alignment##
###########################
"""
Will read in an alignment (-f) in any format (-i),
and rename the strains (-n).  Will output the translated
strain alignment to a new file.
"""

#import modules
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus
import argparse
import os


#parse command line input
PARSER = argparse.ArgumentParser(description=
                                 """
                                 Will read in an alignment (-f) in any format
                                 (-i), and rename the strains (-n).  Will
                                 output the translated strain alignment to a
                                 new file.
                                 """)
PARSER.add_argument('-f', '--filename', help=
                    """
                    Name of input file.
                    """
                    , required=True)
PARSER.add_argument('-n', '--names_file', help=
                    """
                    Names file containing original names (col 1),
                    translated names in col 2.
                    """
                    , required=True)
PARSER.add_argument('-i', '--informat', help=
                    """
                    File format of input. e.g., genbank, fasta, phylip, nexus.
                    """
                    , required=True)
PARSER.add_argument('-o', '--outformat', help=
                    """
                    File format of output. e.g., genbank, fasta, phylip,
                    nexus.
                    """
                    , required=True)

ARGS = PARSER.parse_args()

def alignment_read(infile, informat):
    """
    Read in the alignment, rename the strains according to names file.
    """
    with open(infile, 'r') as input_handle:
        alignment = AlignIO.read(input_handle, informat, alphabet=generic_dna)
        return alignment

def names_read(names):
    """
    Read in the tab-delimited names files, store the name translations
    in a list of tuples.
    """
    name_trans = []
    with open(names, 'r') as input_handle:
        for line in input_handle:
            row = line.rstrip().split(",")
            name_trans.append((row[0], row[1]))
    return name_trans

def translate_record_ids(alignment, names):
    """
    Read the alignment stored in memory, translate the record.ids, output
    record.id-translated alignment.
    """
    alignment_trans = alignment
    for record in alignment:
        for name in names:
#             print record.id.split()
            if name[0]+"_" in record.id+"_":
#                 print name[0], record.id
                record.id = name[1]
    return alignment_trans

def write_alignment(alignment_trans, outformat):
    """
    Read in the translated alignment, write this out to file in any
    format.
    """
    with open(os.path.splitext(ARGS.filename)[0]+"_nametrans."+outformat, "w" \
    ) as output_handle:
        if outformat == "nexus":
            alignment_trans = Nexus.Nexus(alignment_trans.format("nexus"))
            alignment_trans.write_nexus_data(output_handle, interleave=False)
        else:
            AlignIO.write(alignment_trans, output_handle, outformat)
        print '\nAlignment with translated strain names written to "'+\
        output_handle.name+'".'

#Execute the functions above.
INPUT_ALIGNMENT = alignment_read(ARGS.filename, ARGS.informat.replace(".", "")\
.lower())

NAMES_TRANSLATED = names_read(ARGS.names_file)

ALIGNMENT_TRANSLATED = translate_record_ids(INPUT_ALIGNMENT, NAMES_TRANSLATED)

write_alignment(ALIGNMENT_TRANSLATED, ARGS.outformat.replace(".", "").lower())
