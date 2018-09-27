#!/bin/bash 

# run this via the command: 
# qsub -t x -v BEDA=bedfilesa.list -v BEDB=bedfiles2.list -v OUTDIR=outdir run_circular_permutation.sh

#$ -N permute
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=5G 
#$ -j y
#$ -pe sharedmem 4

PROJECTDIR=/home/tballing/bioinfsvice/DSB_model/supplementary_files
unset MODULEPATH
. /etc/profile.d/modules.sh 
module load R

tmpdir=~/scratch/$JOB_NAME/$JOB_ID.$SGE_TASK_ID
mkdir -p $tmpdir

window=0
line=`sed -n "$SGE_TASK_ID"p $BEDA`
set $line
abedid=$1
abed=$2

# Preprocess bedfiles: 
# 1. check if file is gzipped and unzip if it is.
# 2. add "chr" to chromosome names if it's not already there
# 3. Extend the loci in filea by the window size given. 
if file --mime-type "$abed" | grep -q gzip$
then
	gzip -dc $abed \
	| sed '/chr/!s/^/chr/' \
	| awk '{s=$2-win; if (s<0) s=0; print $1"\t"s"\t"$3+win}' win=$window \
	| sort -k1,1 -k2,2n > $tmpdir/a.bed
else
	sed '/chr/!s/^/chr/' $abed \
	| awk '{s=$2-win; if (s<0) s=0; print $1"\t"s"\t"$3+win}' win=$window \
	| sort -k1,1 -k2,2n > $tmpdir/a.bed 
fi 

# Loop through the list of files in the BEDB list
while read bline; do 
set $bline
bbedid=$1
bbed=$2
output=$OUTDIR/$abedid.$bbedid.$window
mkdir -p $output

echo "Working on $bbedid    $bbed..."

if file --mime-type "$bbed" | grep -q gzip$ 
then 	
	gzip -dc $bbed \
	| sed '/chr/!s/^/chr/' | cut -f1-3 | sort -k1,1 -k2,2n > $tmpdir/b.bed
else
	sed '/chr/!s/^/chr/' $bbed | cut -f1-3 | sort -k1,1 -k2,2n > $tmpdir/b.bed
fi 

scripts=$PROJECTDIR/scripts/circperm_scripts
echo "Rscript --vanilla --no-save $scripts/run_circular_permutation.R $tmpdir/a.bed $tmpdir/b.bed $output"
Rscript --vanilla --no-save $scripts/run_circular_permutation.R $tmpdir/a.bed $tmpdir/b.bed $output

# Get simple overlap stats
echo "$scripts/get_overlap_stats.sh $tmpdir/a.bed $tmpdir/b.bed $window $output/overlap_stats.txt"
$scripts/get_overlap_stats.sh $tmpdir/a.bed $tmpdir/b.bed $window $output/overlap_stats.txt

done < $BEDB

# clean up 
# rm -r $tmpdir

