#July 4 2014
#dr.mark.schultz@gmail.com
#github: schultzm
#This file will prune all taxa from the tree, retaining only the tips specified in the tip_labels_to_retain file
#uses Biopython to get the subtree
#to do:
#1. allow multiple trees in input file, or multiple input files
#2. allow it to be used to convert trees without pruning (e.g. by not providing a list of strains to retain).
#probably don't need to have the args.tip_labels_to_retain as str().  Instead, do that at the list building step in the function and convert to strings there
#3. make it so the subtree can be detached
#4. make it so the 'inverse' of the supplied taxa can be used

#First, import the modules
import argparse
from Bio import Phylo

#set up the arguments
parser = argparse.ArgumentParser(description = "This file will prune all taxa from the tree, retaining only the tips specified in the tip_labels_to_retain file (specified with the \'-t\' switch).  Uses Biopython to get the subtree.  If specifying \'nexus\' as output format, delete the second semicolon (\';\') after the tree in the output file (using e.g., TextWrangler) before opening in a tree view program (e.g., Figtree). KNOWN BUGS: An example using a Green-genes tree.  In this case tip names were integers not strings and line breaks were in the tree string: 1. This script worked only after first removing all the line breaks so the tree string was on a single line; 2. Then the tree was then converted from newick to nexus in dendroscope (because tip labels were integers not strings, could also do this conversion with biopython etc.); 3. And on the command line to run this script, it worked when using informat and outformat as 'nexus' (might work with newick and dendroscope but haven't tested); 4. The tree would open in Dendroscope (but not in Figtree - not sure why.  Dendroscope seems to cope better with the integer tip labels).")
parser.add_argument('-i', '--input_tree_file', help = 'Input tree file.', required=True)
parser.add_argument('-f', '--input_tree_format', help = '\'nexus\' or \'newick\'', required = True)
parser.add_argument('-o', '--output_tree_format', help = '\'nexus\' or \'newick\'', required = True)
parser.add_argument('-t', '--tip_labels_to_retain', help = 'File containing tree tip names that you would like to retain.', required = True)

args = parser.parse_args()

#container list to store the name of the output tree file for printing at end of run
output_tree = []

#prune (delete) all tips and associated branches except for the ones in the tip_labels_to_retain file
def filter(input, informat, outformat, tip_labels):
	tree = Phylo.read(input, informat)
# 	print tree
	tip_names_original = [tip.name for tip in tree.get_terminals()]
# 	print tip_names_original
	with open(tip_labels, 'r') as tip_names:
		tip_list_input = [line.rstrip('\n') for line in tip_names]
		tips_to_prune = [i for i in tip_names_original if i not in tip_list_input]
# 		print tips_to_prune
		for i in tips_to_prune:
			tree.prune(i)
		# print str(tree)
		with open(input+'_subtree_pruned.'+outformat+'.tre', 'w') as output_handle:
			output_tree.append(str(output_handle))
			Phylo.draw_ascii(tree)
			Phylo.write(tree, output_handle, outformat)

#execute the script
filter(args.input_tree_file, args.input_tree_format.replace('.',''), args.output_tree_format.replace('.',''), str(args.tip_labels_to_retain))

#Be polite to the user and tell them where their output file went.
print ""
print 'Done'
print 'File written to:'
x = [i.split("'") for i in output_tree]
for i in x:
	print "'"+i[1]+"'"
print ""

