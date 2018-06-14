
.libPaths("/exports/igmm/eddie/NextGenResources/software/R/x86_64-pc-linux-gnu-library/3.2")

library("regioneR") 

args <- commandArgs(trailingOnly=TRUE)
abedfile=args[1]
bbedfile=args[2]
output=args[3]

abed=import.bed(abedfile) 
bbed=import.bed(bbedfile) 

set.seed(123) 
hg19_masked=getBSgenome("BSgenome.Hsapiens.UCSC.hg19.masked") 
abed=filterChromosomes(abed, organism="hg", chr.type="canonical", keep.chr=seqnames(hg19_masked))
bbed=filterChromosomes(bbed, organism="hg", chr.type="canonical", keep.chr=seqnames(hg19_masked))

pt2 <- permTest(A=abed, B=bbed, ntimes=1000, alternative="auto", verbose=FALSE, genome=hg19_masked, evaluate.function=numOverlaps, randomize.function=circularRandomizeRegions, force.parallel=NULL, mc.cores=4, count.once=FALSE, mc.set.seed=FALSE) # mc.se.seed=FALSE needs to be set for reproducing

save(pt2, file=paste(output, "permute_results.Robject", sep="/"))
capture.output(print(pt2), file=paste(output, "permute_summary.txt", sep="/"), type="output")
capture.output(summary(pt2[[1]]), file=paste(output, "circular_results.txt", sep="/"), type="output")


