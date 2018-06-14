#!/bin/bash 

# Run this via the following: 
# qsub -v BEDGRAPH=file.bedgraph.gz -v CHROMSIZES=hg19_chrom.lengths bedgraph_to_bigwig.sh 

#$ -N bgtobw 
#$ -cwd
#$ -j y 
#$ -l h_rt=4:00:00
#$ -l h_vmem=10G
# $ -pe sharedmem 1 

LC_COLLATE=C


tmpdir=$SCRATCH/$JOB_ID
mkdir -p $tmpdir

chromsizes=$CHROMSIZES

f=$BEDGRAPH
label=`basename $f .bedgraph.gz`
zcat -dc $f | sed '/chr/!s/^/chr/' | grep -v ^chrM | sort -k1,1 -k2,2n > $label.sorted.bedGraph
echo "bedGraphToBigWig $label.sorted.bedGraph $chromsizes $label.bigWig"
bedGraphToBigWig $label.sorted.bedGraph $chromsizes $label.bigWig
gzip $label.sorted.bedGraph 


