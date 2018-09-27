#!/bin/bash 

#########################################################
# Set datasets and parameters
PROJECTDIR=/home/tballing/bioinfsvice/DSB_model/supplementary_files #./
scripts=$PROJECTDIR/scripts/circperm_scripts

###########################################
# Launch job to permute the hotspots and
# look for enrichment in different cancer gene sets
infile=$1
infile=`readlink -e $infile`
genesets=$2
genesets=`readlink -e $genesets`

wdir=`basename $infile .sets`
mkdir -p $wdir
cd $wdir
annlist=`basename $genesets .list`
outd=circperm.$annlist
mkdir -p $outd

# make bedfiles from the infile
subsets=`sed 1d $infile | cut -f11 | sort -u`
for s in $subsets; do
	grep $s $infile | awk '{printf "%s\t%i\t%i\n", $1, $2, $3}' > $s.bed
	ls *.bed | awk '{print $1"\t"$1}' | sed 's/.bed//1' > beds.list
done 

n=`wc -l < $genesets`
qsub -t 1-$n -v BEDA=$genesets -v BEDB=beds.list -v OUTDIR=$outd $scripts/run_circular_permutation.sh


