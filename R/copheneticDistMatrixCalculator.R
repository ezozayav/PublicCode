#!/usr/bin/env Rscript
ptm <- proc.time()
help = cat(
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Will find the pairwise distance matrix between tips in a tree file.  
The distance matrix columns and rows will be ordered by tip name.  

Arguments in this order:
    1) treefile.tree        (required, name of the tree file)
    2) treeformat           ('nexus' or 'newick')

Example usage:
#first do make the R script executable by doing:
chmod 777 copheneticDistMatrixCalculator.R
#then execute it by doing:
./copheneticDistMatrixCalculator.R treefile.tree treeformat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\n\n")

#Read options from the command line
args = commandArgs(trailingOnly=TRUE)

#print the help
if(is.element("--help", args) | is.element("-h", args) | is.element("-help", args) | is.element("--h", args)){
    cat(help,sep="\n")
    stop("\nHelp requested.")
}

#ensure that two arguments are used
if(length(args)<2) {
    cat(help,sep="\n")
    stop("\nIncorrect usage, see help above.  Tree file and tree format required.\n")
}

#process command line arguments:
tree_file=args[1]
tree_fmt=args[2]

#print the arguments provided by the user
print(paste("Tree file:", tree_file))
print(paste("Tree format:", tree_fmt))

#stop if invalid tree format given
#if(tree_fmt != "newick" | tree_fmt != "nexus"){
#    stop("\nThe tree file must be exactly 'nexus' or 'newick'.")}

readtree = function(treefile, treeformat){
    if(treeformat=="newick"){
    t = read.tree(treefile)
    return(t)
    }
    if(treeformat=="nexus"){
    t = read.nexus(treefile)
    return(t)
    }
}

library(ape)
tree = readtree(tree_file, tree_fmt)
#to save the tree plot to file uncomment the next three lines
#pdf("test.tree.plot.pdf")
#plot.phylo(tree)
#dev.off()
cophen_dist = cophenetic.phylo(tree)
cophen_dist_ordered=cophen_dist[order(colnames(cophen_dist)),order(colnames(cophen_dist))]
#to preview the distance matrix and see its class, uncomment the next two lines
#head(cophen_dist)
#class(cophen_dist)
out_file = paste(tree_file, ".dist.tab", sep="")
write.table(cophen_dist_ordered, out_file, row.names = TRUE, sep=",", na = "", qmethod="double")
#to read it in the newly created out_file for further analysis do:
#m=read.table(outfile, header=T, sep=",", row.names=1, check.names=F, as.is=T)

cat("Processing completed in (seconds):\n")
proc.time() - ptm
cat("\n")
