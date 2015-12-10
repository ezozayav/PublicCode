#May 16 2014
#updated Sept 18 2014 to input and output any tree format from newick or nexus.
#NB: the nexus output has two semicolons at the end of the tree.  Remove one manually
#before opening in a tree viewer
#dr.mark.schultz@gmail.com
#github: https://github.com/schultzm
#reads in unlimited number of trees.  

#First, import the modules	
import argparse	
from Bio import Phylo	
import glob	
	
#set up the arguments parser to deal with the command line input	
parser = argparse.ArgumentParser(description = "Takes a set of tree files, which will be searched for based on files in the working directory with that extension (-e switch), and then renames the tips based on names in a tab-delimited namesfile.")	
parser.add_argument('-e', '--filesuffix', help = "e.g., '.tree'", required = True)
parser.add_argument('-i', '--treeformat_input', help = "nexus or newick", required = True)
parser.add_argument('-o', '--treeformat_output', help = "nexus or newick", required = True)
parser.add_argument('-n', '--namesfile', help = "tab delimited file with original names in column 1, new names in column 2", required = True)
args = parser.parse_args()	
	
def readconvert(filesuffix, treeformat_input, treeformat_output, namesfile):
	idtable = {}
	f = open(namesfile, "r")
	for line in f:
		fields = line.rstrip().split("\t")
		idtable[fields[0]] = fields[1]
	#this is the list containing the file names
	filelist = glob.glob('*.'+str(filesuffix.replace('.','')))
	for i in filelist:
		tree = Phylo.parse(i, treeformat_input)
		for t in tree:
			for node in t.get_terminals():
				name = node.name
				if name in idtable:
					node.name = idtable[name]
				else:
					node.name = name
					print name +' not in table'
			Phylo.write(t,i.replace('.tree', '_tipsrenamed.tree'), treeformat_output)
readconvert(args.filesuffix, args.treeformat_input, args.treeformat_output, args.namesfile)
