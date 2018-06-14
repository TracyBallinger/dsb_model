#!/bin/sh 

#######################################
# Set variables for directories, etc. 

# ENCODEDIR: where the .bigwig files from encode will be downloaded.
# There should be a lot of space here.
ENCODEDIR=~/scratch
scripts=$PROJDIR/scripts
fscripts=$scripts/feature_scripts

############################################################
# Download UCSC bigWig files
cd $ENCODEDIR
mkdir -p NHEK_UCSC
cd NHEK_UCSC
qsub -v FILES=$fscripts/Nhek_ext_bigwigs.list $fscripts/download_ucscfiles.sh
cd ../ 
mkdir -p K562_UCSC
cd K562_UCSC
qsub -v FILES=$fscripts/K562_bigwigs.list $fscripts/download_ucscfiles.sh 
cd ../
mkdir -p MCF7_UCSC
cd MCF7_UCSC
qsub -v FILES=$fscripts/MCF7_bigwigs.list $fscripts/download_ucscfiles.sh 

