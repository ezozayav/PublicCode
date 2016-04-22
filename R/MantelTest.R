#!/usr/bin/env Rscript
ptm <- proc.time()
help = cat(
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Performs a mantel test between two csv-delimited matrices.  Select the number 
of permutations.  Arguments in this order:
    1) distance matrix 'a' (required, name of matrix 'a')
    2) distance matrix 'b' (required, name of matrix 'b')
    3) correlation method  (required, either 'pearson', 'spearman' or 'kendall',
                            NB: 'kendall' is slow to calculate)
    4) n permuations       (required, '999' is typically used as default)
    
Example usage:
#first, render the R script executable by doing:
chmod 777 MantelTest.R
#then execute it by doing:
./MantelTest.R matrix_a matrix b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\n\n")

#Read options from the command line
args = commandArgs(trailingOnly=TRUE)

#print the help
if(is.element("--help", args) | is.element("-h", args) | is.element("-help", args) | is.element("--h", args)){
    cat(help,sep="\n")
    stop("\nHelp requested.")
}

#ensure that four arguments are used
if(length(args)<4) {
    cat(help,sep="\n")
    stop("\nIncorrect usage, see help above.\n")
}

#process command line arguments:
matrix_a = args[1]
matrix_b = args[2]
corr_method = args[3]
n_perms = as.numeric(args[4])

#print the arguments provided by the user
print(paste("Matrix 'a':", matrix_a))
print(paste("Matrix 'b':", matrix_b))
print(paste("Correlation method:", corr_method))
print(paste("Number of permutations:", n_perms))

library(vegan)

ma = read.table(matrix_a, header=T, sep=",", row.names=1, check.names=F, as.is=T)
mb = read.table(matrix_b, header=T, sep=",", row.names=1, check.names=F, as.is=T)

sink(file=paste("MantelTest_", matrix_a, "_vs_", matrix_b, ".txt", sep = ""))
cat("Comparison\tMethod\tStatistic\tSignificance\tPermutations\n")
m_test = mantel(ma, mb, corr_method, n_perms)
#m_test_names = names(m_test)
#"call", "method", "statistic", "signif", "perm", "permutations" "control" 
cat(paste(paste(matrix_a, "vs", matrix_b, sep=" "), m_test$method, 
      m_test$statistic, m_test$signif, m_test$permutations, sep="\t"))
cat("\n")
sink()

cat("Processing completed in (seconds):\n")
proc.time() - ptm
cat("\n")
