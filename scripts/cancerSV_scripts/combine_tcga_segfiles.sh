#!/bin/bash 

# run this via the command: 
# qsub -v OUT=outdir combine_tcga_segfiles.sh   

#$ -N combine 
#$ -cwd
#$ -l h_rt=02:00:00
#$ -l h_vmem=2G
#$ -j y 

unset MODULEPATH
. /etc/profile.d/modules.sh

wdir=~/bioinfsvice/supplementary_files/data/cancer_SVcnts/tcga
tmpdir=~/scratch/randbkpt
mkdir -p $tmpdir
tmpdir=$tmpdir/$JOB_ID
mkdir -p $tmpdir

outdir=`readlink -e $OUT`
Nextgendir=~/NextGenResources

# Create .seg file for each tissue type and also for amplifications and deletions separately 
cancertypes=(`cut -f2 $wdir/TCGA_SV_calls.info.txt | sort -u`)
for ctype in ${cancertypes[@]}; do 
	echo "working on $ctype...."
	grep $ctype $wdir/TCGA_SV_calls.info.txt > tmplist
	rm tmp.segs.txt
	while read line; do 
		set $line
		f=`grep $3 $wdir/TCGA_SV_files.list`
		if [[ ! -z $f ]]; then 
			cat $f >> tmp.segs.txt 
		else
			echo "$ctype $3 not found" >> unfound_files.list
		fi
	done < tmplist
	head -1 tmp.segs.txt > tmph
	grep -v ^Sample tmp.segs.txt | cat tmph - > $outdir/$ctype.segs.txt
done 

#########################################################
# In these files, every part of the genome seems to have a segment called, 
# most with a mean close to zero.  We just want the segments that
# definitely have a copy number change for amplifications and deletions.  

for ctype in ${cancertypes[@]}; do
	sed 1d $outdir/$ctype.segs.txt | awk '$6 < -1' \
	| cat tmph - > $outdir/$ctype.del.txt
	sed 1d $outdir/$ctype.segs.txt | awk '$6 > 1' \
	| cat tmph - > $outdir/$ctype.amp.txt
done 

