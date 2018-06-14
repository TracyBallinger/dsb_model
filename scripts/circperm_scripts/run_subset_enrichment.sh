#!/bin/bash 

#########################################################
# Set datasets and parameters
PROJECTDIR=./
scripts=$PROJECTDIR/scripts/circperm_scripts

###########################################
# Launch job to permute the hotspots and
# look for enrichment in different cancer gene sets
infile=$1
genesets=$2

wdir=`basename $infile .sets`
mkdir -p $wdir
n=`wc -l < $genesets`
# make bedfiles form the infile
subsets=`sed 1d $infile | cut -f11 | sort -u`
for s in $subsets; do
	cd $wdir
	grep $s ../$infile | awk '{printf "%s\t%i\t%i\n", $1, $2, $3}' > $s.bed
	ls *.bed | awk '{print $1"\t"$1}' | sed 's/.bed//1' > beds.list
	mkdir -p circperm
	qsub -t 1-$n -v BEDA=$genesets -v BEDB=beds.list -v OUTDIR=circperm $scripts/run_circular_permutation.sh
	cd ../
done


