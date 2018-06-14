#!/bin/bash

# qsub -v FILES=filelist 

#$ -cwd
#$ -j y
#$ -N wget
#$ -l h_rt=8:00:00
#$ -l h_vmem=1G


# Download a bunch of data from here: 

while read line; do 
	fn=$line
	if [ ! -f $fn ]; then  
		rsync -avzP rsync://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/$fn .
	fi  
done < $FILES
