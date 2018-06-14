#!/bin/bash 

# Run this via: 
# qsub -v FEATURES=features_matrix.txt -v BEDGRAPH=data.bedgraph -v OUTDIR=outdir run_randomForest.sh 

#$ -N randforest
#$ -cwd 
#$ -j y 
#$ -l h_rt=8:00:00
#$ -l h_vmem=20G

# Need R 
# module load R

mkdir -p $OUTDIR 
Rscript run_randomForest.R $FEATURES $BEDGRAPH $OUTDIR 

