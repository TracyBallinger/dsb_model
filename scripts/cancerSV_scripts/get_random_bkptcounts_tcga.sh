#!/bin/bash 

# run this via the command: 
# qsub -t 1-n -v SEGLIST=segfiles.list -v BINS=hg19_50kb_bins.bed -v OUT=outdir get_random_bkptcounts_tcga.sh 

#$ -N randbkpt 
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G
#$ -j y 
#$ -v R_LIBS_USER=/exports/igmm/eddie/NextGenResources/software/R/x86_64-pc-linux-gnu-library/3.2

PROJECTDIR=~/bioinfsvice/DSB_model/supplementary_files

unset MODULEPATH
. /etc/profile.d/modules.sh
module load R
module load igmm/libs/libpng/1.2.56 

tmpdir=~/scratch/randbkpt
mkdir -p $tmpdir
tmpdir=$tmpdir/$JOB_ID.$SGE_TASK_ID
mkdir -p $tmpdir

ITERS=10

segfile=`sed -n "$SGE_TASK_ID"p $SEGLIST`
name=`basename $segfile .txt`
outfn=$OUT/$name.rcnts
# the get_hotspots.R takes a bed file as input 
# (Note that the second sed command here adds a chr to the beginning of the line if there isn't already one there. 
sed 1d $segfile | cut -f2-5 | sed '/chr/!s/^/chr/' | LC_ALL=C sort -k1,1 -k2,2n > $tmpdir/segs.bed 

# Need to liftover the segments from hg38 to hg19
Nextgendir=~/NextGenResources
$Nextgendir/software/liftOver $tmpdir/segs.bed $Nextgendir/liftOverResources/hg38ToHg19.over.chain.gz $tmpdir/segs.hg19.bed $tmpdir/seg.hg19.unmapped

Rscriptdir=$PROJECTDIR/scripts/cancerSV_scripts 
Rscript --vanilla --no-save $Rscriptdir/get_hotspots_smallsegs.R $tmpdir/segs.hg19.bed $BINS $outfn $ITERS

