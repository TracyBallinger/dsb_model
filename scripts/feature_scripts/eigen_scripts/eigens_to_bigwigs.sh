#!/bin/bash 

# Run this via the following:
# qsub -v EIG=NHEK.eigens.txt $sedir/eigens_to_bigwigs.sh 

#$ -N bgtobw 
#$ -cwd
#$ -j y 
#$ -l h_rt=4:00:00
#$ -l h_vmem=4G
LC_COLLATE=C

# Make a new file for chrom sizes because the eigenvectors are rounded to the nearest 5Mb and extend past the chrom lengths. bedGraphToBigWig doesn't like it if the data extends past the chromosome lengths. 
awk '{x=int($2/50000); newl=50000 * (x+1); print $1"\t"newl}' hg19_chr.lengths \
	| sed 's/^/chr/' | LC_ALL=C sort -k1,1 > chrlengths.txt

label=`basename $EIG .txt`
sort -k1,1 -k2,2n $f > $label.bedGraph
bedGraphToBigWig $label.bedGraph chrlengths.txt $label.bigWig
rm $label.bedGraph

