#Author: Mark Schultz
#email: dr.mark.schultz@gmail.com
#github: schultzm
#library(Rmpi)
#library(snow)
#cl <- makeCluster(8, type = "MPI")
#stopCluster(cl)

directory <- "~/Desktop/HoltBacPathGenomicPostdoc/Acinetobacter/Vietnam/Vietnam_Fig2_GC2SubAlignmentPhylogenies_SNPaware/50_SNPs_s25/MantelTests_50SNPs_s25/"
filename <- "AllTrees_50SNPs_s25.nex.trees_MantelTestResults_correlations_OutgroupExcluded.csv"

setwd(directory)
#read in csv, store it in x, but convert it to a matrix, otherwise the logicals comparison won't work
x <- as.matrix(read.csv(file = filename, header = T, row.names=1))
#edit the row and column names
colnames(x) <- gsub("_matrix", "", colnames(x))
colnames(x) <- gsub("site", "s", colnames(x))
rownames(x) <- gsub("_matrix", "", rownames(x))
rownames(x) <- gsub("site", "s", rownames(x))
#get a matrix of logicals for the lower triangle, with the lower being T
lower_tri_logicals <- lower.tri(x, diag = FALSE)
#compare the logicals to x and store in lower_triangle_x if the logical is T, and store NA if F
lower_triangle_x <- ifelse(lower_tri_logicals==T, x, NA)

#give this new matrix the same row and column names
colnames(lower_triangle_x) <- colnames(x)
rownames(lower_triangle_x) <- rownames(x)

#now grab the values stored at each row and column index, skipping if the value is NA.  
y <- matrix(0, 0, 3)
#the name "Weight" is given to correlation because of the potential to use this output file in Cytoscape.
colnames(y) <- c("Node_1", "Node_2", "Weight")
output <- paste(gsub(".csv", "", filename),"_ColumnisedMatrix.tab.txt", sep = "")
out <- file(output, "w")
cat(colnames(y), file = out, sep="\t")
cat("\n", file =out)

for (i in 1:length(colnames(lower_triangle_x))){
	for (j in 1:length(rownames(lower_triangle_x))){
		if (is.na(lower_triangle_x[i,j])) next
		cat(c(colnames(lower_triangle_x)[i],rownames(lower_triangle_x)[j],(lower_triangle_x[i,j])), file=out, sep="\t")
		cat("\n", file =out)
		}}

close(out)

#now read in the output file, return summary stats, histogram and box-plot of the weights to decide on where to set the cutoff.  
corr_matrix <- read.table(file=output, header = T)
summary(corr_matrix[,"Weight"])
pdf(file=paste(output, ".pdf", sep=""))
hist(corr_matrix[,"Weight"], breaks=length(corr_matrix[,"Weight"])/1000,)
boxplot(corr_matrix[,"Weight"])
dev.off()
