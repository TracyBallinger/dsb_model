#!/bin/bash 

# run this via:
# qsub -v OUTDIR=circdir -v OUTFILE=outfile.dat gather_permute_data.sh  

#$ -N gathercp
#$ -cwd
#$ -l h_vmem=2G
#$ -j y 

OUTDIR=$1
OUTFILE=$2

echo -e "bedA\tbedB\tpvalue\tzscore\talt\tregionset\twindow\tAcnt\tBcnt\tAbp\tBbp\tintersection\tunion\tjaccard\tNinter\tAinter\tBinter\tAuniq\tBuniq\tA_fifty\tB_fifty" > $OUTFILE

for f in `ls $OUTDIR/*/permute_summary.txt`; do 
d=`dirname $f`
namea=`echo $d | sed "s&$OUTDIR/&&" | awk '{split($1, a, "."); print a[1]}'`
nameb=`echo $d | sed "s&$OUTDIR/&&" | awk '{split($1, a, "."); print a[2]}'`

pvalue=`grep ^P-value $d/permute_summary.txt | awk '{print $2}'`
zscore=`grep ^Z-score $d/permute_summary.txt | awk '{print $2}'`
alternative=`grep ^Alternative $d/permute_summary.txt | awk '{print $2}'`
region=`grep region $d/permute_summary.txt | awk '{print $7}'`

stats=`cut -f3- $d/overlap_stats.txt | sed 1d` 
echo -e "$namea\t$nameb\t$pvalue\t$zscore\t$alternative\t$region\t$stats" >> $OUTFILE 
done 

