tree_file <- "CP001921_alleles_VNClone2subset_3outgroups_n25_19509SNPs_var_cons0.90_var.nexus_INCLUDEsnps_INCLUDE_include_SNPs_recategorised_Correct_MarksMethod.csv.nexus.tre_subtree_pruned.newickDATErootedNoHash.tre"
tree_format <- "newick"
library(adephylo)
library(ape)
if(tree_format=="newick")
{
	tree <- read.tree(tree_file)
}
par(mfrow=c(1,2))
plot(tree)
tips <- tree$tip.label
dist_dates <- matrix(,nrow=0,ncol=2)
colnames(dist_dates) <- c("year", "tip_to_root_distance")
dist_dates

for (i in tips)
{
	y <- unname(distRoot(tree, i, "patristic"))
	x <- as.numeric(strsplit(gsub("'","",i), "_")[[1]][5])
	dist_dates <- rbind(dist_dates,c(x,y))

}
linear_model <- lm(dist_dates[,"tip_to_root_distance"]~dist_dates[,"year"])
linear_model_coeff <- round(coef(linear_model), 4)
#linear_model_coeff
#linear_model_coeff[2]
#linear_model_coeff[1]
x_int <- round(linear_model_coeff[1]*-1/linear_model_coeff[2], 4)
r2 <- format(summary(linear_model)$adj.r.squared, digits=4)
corr_coeff <- round(cor(dist_dates[,"tip_to_root_distance"], dist_dates[,"year"]),4)
slope <- round(linear_model_coeff[2], 4)
date_range <- max(dist_dates[,"year"])-min(dist_dates[,"year"])
#date_range
#slope
#r2

#print(linear_model_coeff)
plot(y=dist_dates[,"tip_to_root_distance"], x=dist_dates[,"year"], xlab="Date of isolation", ylab="Root-to-tip divergence (patristic distance)"); abline(linear_model, col="red")
legend("bottomright", bty="n", legend=c(paste("MRCA (x-int) year is ", x_int), paste("R2 is", r2 ), paste("Correlation coefficient is ", corr_coeff), paste("Slope (rate) is ", slope), paste("Date-range is ", date_range, "years")))
mtext(bquote(y == .(linear_model_coeff[2])*x + .(linear_model_coeff[1])), 
adj=1, padj=0)




