#!/bin/bash 

# run this via: 
# qsub -t n -v BEDFILES=bedpairs.list -v WINDOW get_overlap_stats.sh  

#$ -N overlap
#$ -cwd
#$ -l h_vmem=3G
#$ -l h_rt=00:30:00
#$ -j y

function get_bed_stats {
    abed=$1
    bbed=$2
    atot=(`awk 'BEGIN{tot=0; cnt=0}{tot +=($3-$2+1); cnt++}END{print cnt"\t"tot}' $abed`)
    btot=(`awk 'BEGIN{tot=0; cnt=0}{tot +=($3-$2+1); cnt++}END{print cnt"\t"tot}' $bbed`)
    jstats=`bedtools jaccard -a $abed -b $bbed | sed 1d`
    aintcnt=`bedtools intersect -u -a $abed -b $bbed | wc -l`
    bintcnt=`bedtools intersect -u -a $bbed -b $abed | wc -l`
    avcnt=`bedtools intersect -v -a $abed -b $bbed | wc -l`
    bvcnt=`bedtools intersect -v -a $bbed -b $abed | wc -l`
	arovrlp=`bedtools intersect -f 0.50 -r -a $abed -b $bbed | wc -l`
	brovrlp=`bedtools intersect -f 0.50 -r -a $bbed -b $abed | wc -l`
    echo -e "${atot[0]}\t${btot[0]}\t${atot[1]}\t${btot[1]}\t$jstats\t$aintcnt\t$bintcnt\t$avcnt\t$bvcnt\t$arovrlp\t$brovrlp"
}

abed=$1
abedid=`basename $abed .gz`
bbed=$2
bbedid=`basename $bbed .gz`
window=$3
output=$4  

echo -e "A\tB\twindow\tAcnt\tBcnt\tAbp\tBbp\tintersection\tunion\tjaccard\tNinter\tAinter\tBinter\tAuniq\tBuniq\tA_fifty\tB_fifty" > $output

stats=`get_bed_stats $abed $bbed`
echo -e "$abedid\t$bbedid\t$window\t$stats" >> $output
 
