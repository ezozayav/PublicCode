#!/usr/bin/env Rscript
# make this script executable by doing 'chmod +x MaskAllelesCSVsitesClonalML.R'
# dr.mark.schultz@gmail.com, github: schultzm
# date: "11 September 2015"
ptm <- proc.time()
help = cat(
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MaskAllelesCSVsitesClonalML.R takes the recombinant blocks output from the Clo-
nalFrameML analysis and masks those sites in the alleles.csv file. The x_allel-
es.csv must be a comma-separated file with the first column titled 'Pos', which
contains the SNP positions. The recombinant blocks file (i.e., the x_output.im-
portation_status.txt file) is a tab-separated file with the first column titled
'Node' containing node names (exactly as in the x_alleles.csv file when the no-
de is a sample name) and, the second and third columns titled 'Beg' and 'End',
respectively (i.e., the beginning and end of the ClonalFrameML-detected recomb-
inant blocks).  In masking, this script will output a copy of the x_alleles.csv
with sites replaced by characters specified after --arg3 (see below).

    Arguments in this order:
    1) x_alleles.csv        (required, name of the alleles.csv file, comma-delimited)
    2) blocks.tsv           (required, name of the ClonalFrameML 'blocks' file,
                                 tab-delimited)
    3) mask-character       (required, character to mask sites with.  E.g., 'N', '?')
    4) Existing
       gap-character        (optional, e.g., '-', '?','N'.  If not given, existing gaps will not be replaced with the masking character.)

Example usage:
./MaskAllelesCSVsitesClonalML.R x_alleles.csv x.output.importation_status.txt ? Y

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\n\n")

# Read options from command line
args = commandArgs(trailingOnly = TRUE)
if(is.element("--help", args) | is.element("-h", args) | is.element("-help", args) | is.element("--h", args)){
    cat(help,sep="\n")
    stop("\nHelp requested.")
}

if(length(args)<3) {
    cat(help,sep="\n")
    stop("\nIncorrect usage, see help above.  At least three arguments required.\n")
}

# Process the command line arguments:

alleles_file = args[1]
blocks_file = args[2]
maskChar = args[3]
if(nchar(maskChar)> 1) {
    stop("\nThe masking character must be a single character. ")
}

existingGapChar = ifelse(length(args)==4,args[4],NA)
if((nchar(existingGapChar) > 1) & !is.na(existingGapChar)) {
    stop("\nThe existing gap character must be a single character. ")
}
cat("Start log...\n")
cat("Arguments passed from the command line to this job:\n")
cat("'x_alleles.csv': ", alleles_file, "\n")
cat("'blocks.tsv': ", blocks_file, "\n")
cat("Recombinant characters in 'x_alleles.csv' to be masked with '", maskChar, "'\n")
if(!is.na(existingGapChar)){
    cat("Will also attempt to replace existing '", existingGapChar, "' in the masked 'x_alleles.csv' with '", maskChar, "'.\n")
}
cat("----\n")

# Read in and preview the data:
alleles = read.csv(file = alleles_file, header=TRUE, check.names = FALSE, as.is = TRUE, sep = ",")
#turn on following two lines for debugging
#cat("\n", alleles_file, "('x_alleles.csv') looks like this:", "\n")
#head(alleles)
blocks = read.csv(file = blocks_file, header=TRUE, check.names = FALSE, as.is = TRUE, sep = "\t")
blocks[,"Beg"] = as.integer(blocks[,"Beg"]) # Convert the block positions to integers
blocks[,"End"] = as.integer(blocks[,"End"]) # Convert the block positions to integers
#turn on following two lines for debugging
#cat("\n", blocks_file, "('blocks.tsv') looks like this:", "\n")
#head(blocks)

#find which Nodes of the blocks file are isolate names in alleles
blocks=blocks[blocks$Node %in% colnames(alleles),]

#split up the blocks into isolates to reduce loop sizes
blocks=split(blocks, blocks$Node)

#go through each split blocks and replace alleles with pos between columns
#2 and 3 of the block with the maskChar
for(i in 1:length(names(blocks))){
    working_blocks = blocks[i][[1]]
    isolate=working_blocks[1,1]
    for(j in 1:nrow(working_blocks)){
        alleles[which(alleles[,1]>=working_blocks[j,2]&alleles[,1]<=working_blocks[j,3]),isolate]=maskChar
    }
}

#replace user-defined gaps with the maskChar
alleles[alleles==existingGapChar] = maskChar

## write out the data
out_file = paste(gsub(".csv", "", alleles_file), "_masked.csv", sep="")
#cat("Saving output to:", out_file, "\n")
write.csv(alleles, out_file, row.names = FALSE, na = "", quote = FALSE)
cat("Processing completed in (seconds):\n")
proc.time() - ptm
cat("\n")
