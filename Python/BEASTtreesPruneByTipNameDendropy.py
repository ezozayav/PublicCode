#!/usr/bin/env python


'''
Prune BEAST trees by tip name using Dendropy in Python
(https://pythonhosted.org/DendroPy/index.html).  Read BEAST trees into
tree_yielder.  For each tree, retain tips with names matching those on
the command line.  Save the pruned trees to file, adding the suffix to
the file name.  Runs slowly on large tree sets.  Use screen mode in
terminal for long runs.

Example usage:
    time python BEASTtreesPruneByTaxonNameDendropy.py --help

email: dr.mark.schultz@gmail.com
github: schultzm
YYYYMMDD_HHMM: 20160912_1000
'''


import argparse
import dendropy


#set up the arguments
PARSER = argparse.ArgumentParser(description='Prune BEAST trees.')
PARSER.add_argument('-f', '--tree_file', nargs='+', help='Path to treefile.',
                    required=True)
PARSER.add_argument('-t', '--taxa', nargs='+', help='Space delimited list.',
                    required=True)
PARSER.add_argument('-s', '--suffix', help='File keyword (for output suffix).',
                    required=True)
ARGS = PARSER.parse_args()


def main(tree_file):
    '''
    Given a BEAST.trees file, prune the tips in the tree to include only those
    tips listed in ARGS.taxa.  Name the output file to include ARGS.suffix.
    '''
    taxon_namespace = dendropy.TaxonNamespace()
    t_file = open(tree_file, "r")
    tree_yielder = dendropy.Tree.yield_from_files(
        files=[t_file],
        schema="nexus",
        taxon_namespace=taxon_namespace,
        store_tree_weights=True,
        preserve_underscores=True,
        rooting="default-rooted",
        ignore_unrecognized_keyword_arguments=True,
        )
    taxa = ARGS.taxa
    subtrees = []
    for tree in tree_yielder:
        tree.retain_taxa_with_labels(taxa)
        subtrees.append(tree)
    trees = dendropy.TreeList(subtrees)
    trees.write(path=tree_file.replace('.trees', ARGS.suffix+'_pruned.trees'),
                schema='nexus')


if __name__ == '__main__':
    for i in ARGS.tree_file:
        main(i)
