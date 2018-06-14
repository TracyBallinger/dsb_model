#!/bin/sh 

#######################################
# Set variables for directories, etc. 

# ENCODEDIR: where the .bigwig files from encode will be downloaded.
# There should be a lot of space here.
ENCODEDIR=~/scratch
PROJDIR=./ 
scripts=$PROJDIR/scripts
fscripts=$scripts/feature_scripts
# FEATURES Where the final feature sets will be put
OUTDIR=./

############################################################
# Score the bedfile with the bigwigs.
# The genomic segments to score with the bigwigs from encode 
refbed=$PROJDIR/hg19_50kb_bins.bed 
# NHEK dataset
awk '{print edir"/NHEK_UCSC/"$1}' edir=$ENCODEDIR NHEK_bigwigs.list > tmpnhek.list 
qsub -v BEDFILE=$bedfile -v BIGWIGS=tmpnhek.list -v OUT=$OUTDIR/NHEK_bigwigs.scores.tab $scripts/score_bedfile_with_bigwigs.sh 
# NHEK extended dataset
awk '{print edir"/NHEK_UCSC/"$1}' edir=$ENCODEDIR NHEK_ext_bigwigs.list > tmpnhek2.list 
qsub -v BEDFILE=$bedfile -v BIGWIGS=tmpnhek2.list -v OUT=$OUTDIR/NHEK_ext_bigwigs.scores.tab $scripts/score_bedfile_with_bigwigs.sh 
# K562 dataset
awk '{print edir"/K562_UCSC/"$1}' edir=$ENCODEDIR K562_bigwigs.list > tmpk562.list 
qsub -v BEDFILE=$bedfile -v BIGWIGS=tmpk562.list -v OUT=$OUTDIR/K562_bigwigs.scores.tab $scripts/score_bedfile_with_bigwigs.sh 
# MCF7 dataset
awk '{print edir"/MCF7_UCSC/"$1}' edir=$ENCODEDIR MCF7_bigwigs.list > tmpmcf7.list 
qsub -v BEDFILE=$bedfile -v BIGWIGS=tmpmcf7.list -v OUT=$OUTDIR/MCF7_bigwigs.scores.tab $scripts/score_bedfile_with_bigwigs.sh 

############################################################
# Clean up 
rm -r $ENCODEDIR
rm tmp*.list
 
 
