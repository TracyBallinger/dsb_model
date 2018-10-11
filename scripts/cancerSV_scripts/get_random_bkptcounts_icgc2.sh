#!/bin/bash 

# run this via the command: 
# qsub -t 1-n -v BEDLIST=bedfiles.list -v BINS=hg19_50kb_bins.bed -v OUT=outdir get_random_bkptcounts_icgc.sh 

#$ -N randbkpt
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=5G
#$ -j y 
#$ -v R_LIBS_USER=/exports/igmm/eddie/NextGenResources/software/R/x86_64-pc-linux-gnu-library/3.2

PROJECTCIR=~/bioinfsvice/DSB_model/supplementary_files

unset MODULEPATH
. /etc/profile.d/modules.sh
module load R

tmpdir=~/scratch/randbkpt
mkdir -p $tmpdir
tmpdir=$tmpdir/$JOB_ID.$SGE_TASK_ID
mkdir -p $tmpdir

ITERS=10
# the get_hotspots.R takes a bed file as input 
# (Note that the second sed command here adds a chr to the beginning of the line if there isn't already one there. 
BED=`sed -n "$SGE_TASK_ID"p $BEDLIST`
outname=`basename $BED .bed`
output=$OUT/$outname.rcnts
cat $BED | sed '/chr/!s/^/chr/' | LC_ALL=C sort -k1,1 -k2,2n > $tmpdir/segs.bed 

Rscriptdir=$PROJECTDIR/scripts/cancerSV_scripts
Rscript --vanilla --no-save $Rscriptdir/get_hotspots_smallsegs.R $tmpdir/segs.bed $BINS $output $ITERS


