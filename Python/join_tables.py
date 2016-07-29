#!/usr/bin/env python

"""
Run CPE_ongoing Nullarbor.pl requests on repeatedly updated list of MDU IDs.
Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20160726_2223

"""

import os
import argparse
import sys
import string
import random
import pymp
from subprocess import call, Popen, PIPE
import shlex
import time
from collections import defaultdict
import tempfile
import shutil
import itertools
from Bio import Phylo
from Bio.Phylo import BaseTree, _io
# from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import _Matrix
from Bio.Phylo.TreeConstruction import _DistanceMatrix
# import Bio.Phylo.TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.PhyloXML import Phylogeny
from Bio.Nexus import Trees
import pandas as pd
# from pandas import *
import networkx, pylab
# from rpy2.robjects.packages import importr
# import rpy2.robjects as ro
# import pandas.rpy.common as com
# from ete2 import PhyloTree #, ProfileFace
# import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import scipy.io as sio

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree, to_tree, leaves_list
import numpy as np
import json
from rpy2.robjects import pandas2ri
# pandas2ri.activate()
# import pandas.rpy.common as com
# import pandas as pd

# from ete2 import Tree
from ete3 import Tree
# import skbio
from multiprocessing import Pool



# set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="CPE ongoing")
PARSER.add_argument("-i", "--mdu_read_ids", help="One MDU-ID per line.",
                    required=True)
PARSER.add_argument("-w", "--wgs_qc", help="Path to WGS\
                    QC. Default '/mnt/seq/MDU/QC/'", 
                    default="/mnt/seq/MDU/QC/", required=False)
PARSER.add_argument("-d", "--delete_tempdirs", help="Delete tempdirs created\
                    during run? Default = 'no'.", default='no', required=False)
PARSER.add_argument("-t", "--threads", help='Number of threads, default=\'72\'',
                    default='72', required=False)
PARSER.add_argument("-a", "--andi_run", help='Run andi phylogenomic analysis?\
                    Default=\'yes\'', default='yes', required=False)
PARSER.add_argument("-m", "--model_andi_distance", help='Substitution model.\
                    \'Raw\', \'JC\', or \'Kimura\'. Default = \'JC\'.', 
                    default='JC', required=False)
PARSER.add_argument('-c', '--percent_cutoff', help='For abricate, call the\
                    gene \'present\' if greater than this value and \'maybe\'\
                    if less than this value. Default = 95. NB: 100 percent = \
                    \'100\', not \'1\'.', default=95, type=int, required=False)

ARGS = PARSER.parse_args()

mlst_schemes = {"Acinetobacter baumannii": "abaumannii",
"Citrobacter freundii": "cfreundii",
"Cronobacter": "cronobacter",
"Enterobacter cloacae": "ecloacae",
"Escherichia coli": "ecoli",
"Klebsiella oxytoca": "koxytoca",
"Klebsiella pneumoniae": "kpneumoniae",
"Pseudomonas aeruginosa": "paeruginosa"}

print mlst_schemes


abricate_headers = []


class Isolate(object):
    def __init__(self, id):
        '''
        Initialise the object with an MDU-ID and a path to the QC files.
        '''
        self.id = id
        self.qc_path = ARGS.wgs_qc
    def assembly(self):
        '''
        Store the path to the assembly, or tell the user if it doesn't exist.
        '''
        isolate_qc_contigs = self.qc_path+self.id+'/contigs.fa'
        if os.path.exists(isolate_qc_contigs):
            return isolate_qc_contigs
        else:
            print isolate_qc_contigs+' does not exist'
    def reads(self):
        '''
        Store the path to the QCd (trimmomatic-ed and FLASHed) reads.
        '''
        pass
    def shortened_id(self):
        '''
        Create a random 9 character (alphanumeric) tag for use as a tempfile
        name.
        '''
        chars=string.ascii_uppercase + string.digits
        return ''.join(random.SystemRandom().choice(chars) for _ in range(10))

    def abricate_path(self):
        abricate_path = self.qc_path+self.id+'/abricate.tab'
        return abricate_path

    def abricate(self):
        '''
        Store the path to abricate results.
        '''
        abricate_path = self.qc_path+self.id+'/abricate.tab'
        ab_data = pd.read_table(abricate_path, sep='\t', header=0)
        nrows = ab_data.shape[0]
        genes = ab_data['GENE'].tolist()
        cov = ab_data['%COVERAGE'].tolist()
        yes = []
        maybe = []
        for i in range(0, nrows):
            if cov[i] >= ARGS.percent_cutoff:
               yes.append(genes[i])
            else:
                maybe.append(genes[i])
        y = {key:'yes' for (key) in yes}
        m = {key: 'maybe' for (key) in maybe}
        ab_results = dict(itertools.chain(y.iteritems(), m.iteritems()))
        abricate_results = pd.DataFrame([ab_results], index=[self.id])
#         print abricate_results
        return abricate_results
        
    def mlst(self, species, assembly):
        '''
        Store the MLST.
        '''
        if species in mlst_schemes:
            cmd = 'mlst --scheme '+mlst_schemes[species]+' --quiet '+assembly
        else:
            cmd = 'mlst --quiet '+assembly
        args_mlst = shlex.split(cmd)
        #should write a proc function to do this call as it's used often
        proc = Popen(args_mlst, stdout = PIPE)
        output = proc.stdout.read()
        #should write a function here too
        if species in mlst_schemes:
            mlst = output.rstrip().split('\n')
            mlst = [line.strip().split('\t') for line in filter(None, mlst)]
            header = mlst[0]
            data = mlst[1]
            ncol = len(header)
            mlst_formatted_dict = {'0_1_scheme': data[1], '0_ST': data[2]}
            k = 1
            for i in range(3, ncol):
                locus = header[i]+'('+data[i]+')'
                mlst_formatted_dict[str(k)] = locus
                k += 1
        else:
            out = output.rstrip().split('\t')[1:]
            ncol = len(out)
            mlst_formatted_dict = {'0_1_scheme': out[0], '0_ST': out[1]}
            k = 1
            for i in range(3, ncol):
                mlst_formatted_dict[str(k)] = out[i]
                k += 1
        mlst_results = pd.DataFrame([mlst_formatted_dict], index = [self.id])
#         print mlst_results
        return mlst_results

    def kraken_reads(self):
        '''
        Get the kraken best hit from reads.
        '''
        #Pipe these commands together
        cmd_grep = "grep -P '\tS\t' /mnt/seq/MDU/QC/"+self.id+"/kraken.tab"
        cmd_sort = 'sort -k 1 -r'
        cmd_head = 'head -2'

        #Split the cmds using shlex, store in args
        args_grep = shlex.split(cmd_grep)
        args_sort = shlex.split(cmd_sort)
        args_head = shlex.split(cmd_head)
        
        #Pipe the output of one args to another
        proc1 = Popen(args_grep, stdout = PIPE)
        proc2 = Popen(args_sort, stdin=proc1.stdout, stdout=PIPE)
        proc3 = Popen(args_head, stdin=proc2.stdout, stdout=PIPE)
        
        output = proc3.stdout.read()
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line in filter(None, kraken)]
        return kraken

    def kraken_contigs(self):
        '''
        Get the kraken best hit from assemblies.
        '''
        #Pipe these commands together
#         kraken --fasta-input /mnt/seq/MDU/QC/2016-15685/contigs.fa | kraken-report | grep -P '\tS\t' | sort -k 1 -r | head -3
        cmd_kraken = 'nice kraken --threads 4 --db /bio/db/kraken/minikraken --fasta-input /mnt/seq/MDU/QC/'+self.id+'/contigs.fa'
        cmd_krk_r = 'kraken-report'
        cmd_grep = "grep -P '\tS\t'"
        cmd_sort = 'sort -k 1 -r'
        cmd_head = 'head -3'

        #Split the cmds using shlex, store in args
        args_kraken = shlex.split(cmd_kraken)
        args_krk_report = shlex.split(cmd_krk_r)
        args_grep = shlex.split(cmd_grep)
        args_sort = shlex.split(cmd_sort)
        args_head = shlex.split(cmd_head)
        
        #Pipe the output of one args to another
        proc1 = Popen(args_kraken, stdout = PIPE)
        proc2 = Popen(args_krk_report, stdin=proc1.stdout, stdout=PIPE, stderr=PIPE)
        proc3 = Popen(args_grep, stdin=proc2.stdout, stdout=PIPE, stderr=PIPE)
        proc4 = Popen(args_sort, stdin=proc3.stdout, stdout=PIPE, stderr=PIPE)
        proc5 = Popen(args_head, stdin=proc4.stdout, stdout=PIPE, stderr=PIPE)
        
        output = proc5.stdout.read()
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line in filter(None, kraken)]
        return kraken

    def species(self, kraken):
        '''
        Based on the results of the kraken scans, define the species.
        '''
        #If the #1 match is more than twice the #2 match, keep #1, else 'indet'
        if len(kraken) > 1:
            if float(kraken[0][0])/float(kraken[1][0]) >= 2:
                species = kraken[0][5].lstrip()
            else:
                species = 'indet'
            return species
        if len(kraken) == 1:
            if float(kraken[0][0]) > 5:
                species = kraken[0][5].lstrip()
            else:
                species = 'indet'
            return species
        
    def species_conclusion(self, kraken_reads_species, kraken_contigs_species):
        krs = kraken_reads_species
        kcs = kraken_contigs_species
        if krs != kcs and krs != 'indet' and kcs != 'indet':
            species_final = 'mixed'
        if krs == 'indet' or kcs == 'indet':
            species_final = 'indet'
        if krs == kcs:
            species_final = krs
        return species_final

    def kraken_percentage(self, kraken):
        '''
        What percentage matched the best hit? Look in the first item of kraken.
        '''
        return kraken[0][0]

def iso_sub_list(id):
    '''
    Takes an id, looks for it in the mdu-wgs.tab file, returns a list of hits.
    '''
    args = shlex.split('grep '+id+' /mnt/seq/MDU/QC/mdu-wgs.tab')
    args_cut = shlex.split('cut -f 1')
    proc1 = Popen(args, stdout = PIPE)
    proc2 = Popen(args_cut, stdin = proc1.stdout, stdout=PIPE)
    output = filter(None, proc2.stdout.read().rstrip().split('\n'))
    return output

def isolate_request_ids(id_file):
    '''
    Reads in the MDU IDs from the request list.
    '''
    request_ids = filter(None, [line.strip() for line in open(id_file).readlines()])
    return request_ids

def make_tempdir():
    '''
    Make a temporary directory.
    '''
    tempdir = tempfile.mkdtemp(dir='.')
    return tempdir

def read_file_lines(file_name):
    '''
    Read in lines of a file, store in a list.
    '''
    file_contents = [line.rstrip().split() for line in open(file_name).readlines()]
    return file_contents
    
def lower_tri(full_matrix):
    '''
    Take a symmetrical matrix, convert it to a lower triangle _Matrix object.
    '''
    lower_tri = []
    names = []
    k = 2
    for i in full_matrix:
        lower_tri.append(map(float, i[1:k]))
        names.append(i[0])
        k += 1
    matrix = _DistanceMatrix(names, lower_tri)
    return matrix

def isolates_available(ids):
    '''
    Of the request IDs, find which ones are actually available for analysis.
    '''
    avail = []
    not_avail = []
    for id in ids:
        isos = iso_sub_list(id)
        if len(isos) == 0:
            not_avail.append(id)
        else:
            for j in isos:
                avail.append(j)
    print '\nOf the isolate IDs requested'
    print 'Found:\n'+'\n'.join(avail)
    print '\nNot found:\n'+'\n'.join(not_avail)
    return avail

def kraken_contigs_multiprocessing(iso):
    id = Isolate(iso)
    kraken_contigs_hit = id.kraken_contigs()
    species_contigs = id.species(kraken_contigs_hit)
    percent_match = id.kraken_percentage(kraken_contigs_hit)
    return iso, species_contigs, percent_match

def main():
    '''
    Read in the MDU-IDs from file. For each id, instantiate an object of 
    class Isolate.  This class associates all of the QC data with the ID tag.  
    Move the contigs for all isolates into a tempdir, with a 9-character 
    filename.  Run andi phylogenomics on all the contig sets.  Infer an NJ tree
    using Bio Phylo from the andi-calculated distance matrix.  Correct the 
    negative branch lengths in the NJ tree using ETE2.  Export the tree to file.
    Export the metadata to file.  Todo: plot the metadata next to tree.'''
    
    ids = isolate_request_ids(ARGS.mdu_read_ids)
    isos = isolates_available(ids)
    
    #make tempdir and store the contigs there for 'andi' analysis
    assembly_tempdir = make_tempdir()
    print '\nMoving contigs into temp place:'
    #Store translation dict (random 9-character filename: original filename
    iso_id_trans = {}
    for iso in isos:
        id = Isolate(iso)
        assembly_path = id.assembly()
        short_id = id.shortened_id()
        iso_id_trans[iso] = short_id
        cmd = 'cp '+assembly_path+' '+assembly_tempdir+'/'+short_id+'_contigs.fa'
        os.system(cmd)
        print cmd
    print '\nTranslated isolate IDs:\nShort\tOriginal'
    for key, value in iso_id_trans.items():
        print value+'\t'+key

    if 'y' in ARGS.andi_run.lower():
        metadata_results = {key: [] for key in isos}
        print 'metadata results:', metadata_results
        base = ARGS.mdu_read_ids.replace('.txt', '')
        andi_mat = 'andi_'+base+'.mat'
        #Run andi
        andi_c = 'nice andi -j -m '+ARGS.model_andi_distance+' -t '+\
                  ARGS.threads+' '+assembly_tempdir+'/*_contigs.fa > '+andi_mat
        print '\nRunning andi with: \''+andi_c+'\''
        os.system(andi_c)
        
        #Read in the andi dist matrix, convert to lower triangle
        dm = read_file_lines(andi_mat)[1:]
        dm = lower_tri(dm)
        
        #Dict to store kraken_contig results
        kraken_contigs_dict = {key: [] for (key) in isos}
        
        #Kraken set at 4 threads, so 18 processes can run on 72 CPUs
        n_isos = len(isos)
        print '\nRunning kraken on the assemblies:'
        if n_isos <= 18:
            p = Pool(n_isos)
            results = p.map(kraken_contigs_multiprocessing, isos)
            for result in results:
                kraken_contigs_dict[result[0]].extend([result[1], result[2]])
        if n_isos > 18:
            p = Pool(18) 
            results = p.map(kraken_contigs_multiprocessing, isos)
            for result in results:
                kraken_contigs_dict[result[0]].extend([result[1], result[2]])

#         print kraken_contigs_dict
        summary_frames = []
        #These step replaces the dm-matrix short names with original names
        for iso in isos:
            iso_df = []
            id = Isolate(iso)
            short_id = iso_id_trans[iso]
            kraken_reads_hit = id.kraken_reads()
            species_reads = id.species(kraken_reads_hit)
            reads_pc = id.kraken_percentage(kraken_reads_hit)
            kraken_contigs_species = kraken_contigs_dict[iso][0]
#             print kraken_contigs_species
            contigs_pc = kraken_contigs_dict[iso][1]
            species_conclsn = id.species_conclusion(species_reads, kraken_contigs_species)
#             print species_reads, reads_pc, kraken_contigs_species, contigs_pc, species_conclsn
            mlst_df = id.mlst(species_conclsn, id.assembly())
            iso_df.append(mlst_df)
#            print 'mlst', mlst_df
            resistome_df = id.abricate()
            iso_df.append(resistome_df)
#            print resistome_df
            kraken_summary_dict = {
                              "__Species_consensus": species_conclsn,
                              "_kraken_reads": species_reads,
                              "_kraken_reads_pc": reads_pc,
                              "_Kraken_contigs": kraken_contigs_species,
                              "_kraken_contigs_pc": contigs_pc
                             }
            kraken_summary_df = pd.DataFrame([kraken_summary_dict], index = [iso])
            iso_df.append(kraken_summary_df)
            metadata_iso = pd.concat(iso_df, axis=1)
            print metadata_iso
            summary_frames.append(metadata_iso) 
            
            for i in range(0, len(dm.names)):
                if dm.names[i] == short_id:
                    #change next line to just 'iso' after debug
                    dm.names[i] = iso#+'__'+species_reads.replace(' ', '_')+'__'+str(reads_pc)
        metadata_all = pd.concat(summary_frames)
        print metadata_all
#         print '\nandi-calculated distance matrix'
#         print dm
#         for i in summary_frames:
#         print help(pd.concat)
#         metadata_results = pd.concat(summary_frames[0:], axis = 0)
#         print metadata_results
        #From the distance matrix, infer the NJ tree
        constructor = DistanceTreeConstructor()
        njtree = constructor.nj(dm)
        njtree.rooted = True
        #Todo: store the temp.tre in a real temp file, not just one called temp.
        Phylo.write(njtree, 'temp.tre', 'newick')
        t = Tree("temp.tre", format=1)
        #Get rid of negative branch lengths
        for node in t.traverse():
            node.dist = abs(node.dist)
        t.set_outgroup(t.get_midpoint_outgroup())
        t.write(format=1, outfile = 'andi_'+base+'BLcorrected.tre')
        print t
        
        ###Development###############
#         dend_fig = 'scipy-dendrogram_'+base+'.pdf'
#         with PdfPages(dend_fig) as pdf:
#             Z = linkage(dm, 'complete')
#             h = dendrogram(Z, color_threshold = 0.6, labels = dm.names, leaf_font_size=2, orientation = 'left') #labels=dm.names
#             clusters = fcluster(Z, 0.6, 'distance')
#             print '\nIsolate\tKrakenSp\tCluster'
#             members_by_cluster = [list(i) for i in zip(dm.names, clusters)]
#             for member in members_by_cluster:
#                 member = [member[0].split('__')[0], member[0].split('__')[1], member[1]]
#             print member
# #         plt.ioff()
# 
#             plt.tight_layout()
#     #         fig = plt.figure()
#             plt.title('Hierarchical Clustering Dendrogram')
#             plt.xlabel('IsolateID_KrakenSp')
#             plt.ylabel('Distance')
#             
#             pdf.savefig()
#             print '\nDendrogram written to '+dend_fig


    if 'y' in ARGS.delete_tempdirs.lower():
        shutil.rmtree(assembly_tempdir, ignore_errors=True)
        print '\nDeleted tempdir '+assembly_tempdir+'\nFinished.\n'

if __name__ == '__main__':
    main()
