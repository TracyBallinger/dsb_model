######################################################
# To run: 
# Rscript run_randomForest2.R features_matrix.txt breakscores.txt outputdir 
# features_matrix.txt and breakscores.txt should have the same number of lines
# features_matrix.txt should have these columns (bed format): 
# <chr><start><end><id><score><scores...>
# breakscores.txt should have these columns (bedgraph format)
# <chr><start><end><score>
# This will build a randomForest using the features to predict the first score
# of the breakscores.txt file


###########################################################
# Read in command line arguments
args <- commandArgs(trailingOnly =TRUE)
featureScoreFile = args[1]
breakscores = args[2]
outputdir = args[3]

###########################################################
# Load appropriate libraries for doing this analysis and plotting
# Need to add this to R's path 
.libPaths("/exports/igmm/eddie/NextGenResources/software/R/x86_64-pc-linux-gnu-library/3.2")
require("randomForest")
library("ggplot2")
# library("gridExtra")
# library("plotrix")
# library("reshape2")
# library("scales")

# setwd("/Users/tballing/Documents/Research/DSB_model")
# all.dat <- read.table("HeLa_bigwig_scores.tab", header=TRUE)
all.dat <- read.table(featureScoreFile, header=TRUE)
# change an NA values to zero 
all.dat[is.na(all.dat)] <- 0 

print(dim(all.dat))
features=all.dat[,5:ncol(all.dat)]

# breakdat=read.table("breakome_scores.tab", header=TRUE)
breakdat=read.table(breakscores, header=TRUE)
features$varx = breakdat[,4]

all.mod <- randomForest(varx ~ ., data=features, ntree=500, do.trace=25, keep.inbag=T, keep.forest=T, nPerm=5, importance=T)
save(all.mod, file=paste(outputdir, "all.mod.rdf", sep="/"))
##############################################################
# Plot which features are important
imp <- importance(all.mod, type=1)
imp <- imp[order(imp[,1], decreasing=T),]
pdf(file=paste(outputdir, "importance.pdf", sep="/"), width=7, height=8)
par(mfrow=c(1,1), mar=c(3,9,1.5,0.5), oma=c(2,0,0,0), mgp=c(0,.5,0))
barplot(rev(imp[1:20]), horiz=T, las=1, cex.names=1.1, col="#FFA50092", border=NA)
mtext(1, outer=T, text="Variable importance (% increase in MSE when permuted)")
dev.off()
write.table(imp, file=paste(outputdir, "variable_importance.txt", sep="/"), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
 
##############################################################
# Plot the profile of the actual values of your variable and the predicted values. 
breaks <- data.frame(cbind(all.dat$chr, all.dat$start, all.dat$end, features$varx, all.mod$predicted))
colnames(breaks) <- c("chr", "start", "end", "breaks", "predicted")
write.table(breaks, file=paste(outputdir, "predicted.txt", sep="/"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")



