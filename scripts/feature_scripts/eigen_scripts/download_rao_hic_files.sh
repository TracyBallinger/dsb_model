#!/bin/bash

#$ -cwd
#$ -j y
#$ -N wget
#$ -l h_rt=4:00:00
#$ -l h_vmem=1G


# Download a hic data from here: 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_NHEK_combined.hic.gz


