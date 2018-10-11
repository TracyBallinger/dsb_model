.libPaths("/exports/igmm/eddie/NextGenResources/software/R/x86_64-pc-linux-gnu-library/3.2")

library(regioneR) 

args <- commandArgs(trailingOnly=TRUE)
segfile = args[1]  # CNV_comparison/Tumorscape/tumorscape_100217.seg.hg19.bed
bbedfile = args[2]
output = args[3]
numiter = as.integer(args[4])

hg19_masked=getBSgenome("BSgenome.Hsapiens.UCSC.hg19.masked")

segs=read.table(segfile, header=FALSE) 
scna=toGRanges(A=segs, format="BED") 
scnafilt <- filterChromosomes(scna, organism="hg19", chr.type="canonical", keep.chr=seqnames(hg19_masked)) 

bins=toGRanges(read.table(bbedfile))

# In doing the counting, if a segment is contained within a genomic bin, I only
# want to count one break, not two so that I'm not double counting, in a sense
# the breakpoints for short segments per bin. 
countOverlap_brks <- function(segs, bins, genome) {
	# get segments that are within a bin. 
	dbcnt <- overlapsAny(toGRanges(segs), bins, type="within")
	nodbls <- segs[!dbcnt,]
	left = data.frame(chr=nodbls[,1], start=nodbls[,2], end=nodbls[,2])
	right = data.frame(chr=nodbls[,1], start=nodbls[,3], end=nodbls[,3])
	dbl <- data.frame(chr=segs[dbcnt,1], start=segs[dbcnt,2], end=segs[dbcnt,3])
	scnabrks=toGRanges(A=data.frame(rbind(left, right, dbl)))
	cnt=countOverlaps(bins, scnabrks)
	return(cnt)
}

randomizeSCNAandCountBrks<- function(scna, bins, genome){
	#perm= circularRandomizeRegions(scna, genome=genome)
	perm=randomizeRegions(A=scna, genome=genome, mask=NA, allow.overlaps=TRUE, per.chromosome=TRUE)
	segs=toDataframe(perm)
	cnt = countOverlap_brks(segs, bins, genome)
	return(cnt)
}

mycnts=matrix(nrow=length(bins), ncol=numiter)
for (i in seq(1:numiter)){
	if ((i %% 50) == 0){
		write(paste("working on iteration", i), stderr())
	}
	mycnts[,i]=randomizeSCNAandCountBrks(scna=scnafilt, bins=bins, genome="hg19")
}

# get the fraction of times that the number of random breakpoints is greater than or equal to the actual number of breakpoints
#segs=read.table(segfile, header=FALSE)
realcnt=countOverlap_brks(segs, bins)
avepcnt <- apply(mycnts, 1, mean)
randomFraction <- apply(cbind(realcnt, mycnts), 1, function(x) sum(x[2:length(x)] > x[1]))/numiter
outdf <- cbind(toDataframe(bins), realcnt, avepcnt, randomFraction)
write.table(outdf, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
# write.table(mycnts, file=paste(output, ".dat", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
