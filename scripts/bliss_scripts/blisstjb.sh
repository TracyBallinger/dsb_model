#!/usr/bin/env bash

# THIS SCRIPT CAN BE CALLED AS
# qsub -v SRR_ID=srr_id -v OUT=outdir blisstjb.sh 
# This pipeline was designed by Silvano Garnerone and modified by Tracy Ballinger. 
# You will need the following other software: 
# scan_for_matches
# bwa
# samtools
# bedtools 

#$ -N bliss
#$ -cwd
#$ -j y 
#$ -l h_rt=02:30:00
#$ -l h_vmem=4G
#$ -pe sharedmem 3

################################################################################
clear
# DEFINING VARIABLES
experiment=$SRR_ID	# e.i. SRR*_1.fastq corresponding to the fastq file name
fastqdir=fastqdir
refgen=GRCh37.fa  	# Reference genome fasta file with bwa index
patfile=example.txt  # this is the pattern file: ../pattern/example.txt
quality=30	 		# mapping quality - they use 30 in the published BLISS paper
numproc=3			# number of threads to be used
scandir=scan_for_matches  # PATH to scan_for_matches

################################################################################
# PREPARE DIRECTORY STRUCTURE
mkdir -p $OUT
mkdir -p $OUT/$experiment
outdir=$OUT/$experiment
tmpdir=$SCRATCH
mkdir -p $tmpdir
tmpdir=$tmpdir/$experiment
mkdir -p $tmpdir

################################################################################
# FIND DATA FILES:
find $fastqdir -type f -iname "*$experiment*.fastq" | sort > filelist_"$experiment"
# PRINT TO TERMINAL THE NAMES OF THE FASTQ FILES THAT HAVE BEEN FOUND
numb_of_files=`cat filelist_"$experiment" | wc -l`
r1=`cat filelist_"$experiment" | head -n1`
echo "R1 is " $r1
if [ $numb_of_files == 2 ]; then
    r2=`cat filelist_"$experiment" | tail -n1`
    echo "R2 is " $r2
fi
rm filelist_"$experiment"

#################################################################
# Reformat files and filter for sequence pattern 
if [ ! -e $tmpdir/r1.2b.aln.fq ]; then 
cat $r1 | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $tmpdir/r1.fa
cat $r1 | paste - - - - | LC_ALL=C sort --temporary-directory=$tmpdir -k1,1 > $tmpdir/r1oneline.fq
# cat $tmpdir/r1.fa | paste - - > $tmpdir/r1oneline.fa
# Filter for the pattern  
cat $tmpdir/r1.fa \
	| parallel --tmpdir $tmpdir --block 100M -k --pipe -L 2 \
	"$scandir/scan_for_matches $patfile - " > $tmpdir/filtered.r1.fa
cat $tmpdir/filtered.r1.fa \
	| cut -d':' -f1 | tr '>' '@' | paste - - \
	| awk '{print $1,$NF}' \
	| LC_ALL=C sort --temporary-directory=$tmpdir -k1,1 > $tmpdir/ID_genomic
LC_ALL=C join $tmpdir/ID_genomic $tmpdir/r1oneline.fq \
	| awk '{print $1,$2,$6,substr($9, length($9)-length($2)+1, length($5))}' \
	| tr " " "\n" > $tmpdir/r1.2b.aln.fq

###################################################################
# If it's paired end, then get the matching second read for the filtered ones.  
if [ ! -z $r2 ]; then
	cat $r2 | paste - - - - \
	| LC_ALL=C sort --temporary-directory=$tmpdir -k1,1 > $tmpdir/r2oneline.fq
	# Filter for the read1s with the pattern  
	LC_ALL=C join $tmpdir/ID_genomic $tmpdir/r2oneline.fq \
	| cut -d' ' -f 1,5,6,9 | tr " " "\n" > $tmpdir/r2.2b.aln.fq 
fi 
fi 

###############################################################
# Align the reads to the reference. 
if [ ! -e $tmpdir/$experiment.q30.sorted.bam ]; then 
if [ -z $r2 ]; then 
	echo "bwa mem -t $numproc $refgen $tmpdir/r1.2b.aln.fq > $tmpdir/$experiment.sam" 
	bwa mem -t $numproc $refgen $tmpdir/r1.2b.aln.fq > $tmpdir/$experiment.sam 
else 
	echo "bwa mem -t $numproc $refgen $tmpdir/r1.2b.aln.fq $tmpdir/r2.2b.aln.fq > $tmpdir/$experiment.sam" 
	bwa mem -t $numproc $refgen $tmpdir/r1.2b.aln.fq $tmpdir/r2.2b.aln.fq > $tmpdir/$experiment.sam 
fi 

###############################################################
# filter out reads with low quality mapping (quality less than 30).
samtools view -Sb -q 30 $tmpdir/$experiment.sam > $tmpdir/$experiment.q30.bam
samtools sort $tmpdir/$experiment.q30.bam -o $tmpdir/$experiment.q30.sorted.bam
samtools index $tmpdir/$experiment.q30.sorted.bam 
# mv $tmpdir/$experiment.q30.sorted.bam $outdir
# mv $tmpdir/$experiment.q30.sorted.bam.bai $outdir
fi 

###############################################################
# Need to get the UMIs for the aligned reads to filter out 
# PCR duplicates downstream 
samtools view -F 0x10 $tmpdir/$experiment.q30.bam \
	| cut -f1,3,4 \
	| awk '{print $0, "\t+"}' > $tmpdir/forward  
samtools view -f 0x10 $tmpdir/$experiment.q30.bam \
	| cut -f1,3,4 \
	| awk '{print $0, "\t-"}' > $tmpdir/reverse
cat $tmpdir/forward $tmpdir/reverse \
	| LC_ALL=C sort --parallel=$numproc -T $tmpdir -k1,1 > $tmpdir/id.chr.loc.strand
cat $tmpdir/filtered.r1.fa | tr -d ">" | cut -d':' -f 1 | paste - - \
	| awk '{if (NF==4) print; else if (NF==5) print $1,$3,$4,$5}' \
	| LC_ALL=C sort --parallel=$numproc -T $tmpdir -k1,1 > $tmpdir/id.umi.barcode.genomic 

LC_ALL=C join $tmpdir/id.chr.loc.strand $tmpdir/id.umi.barcode.genomic \
	| cut -d' ' -f-5 \
	| LC_ALL=C sort --parallel=$numproc -T $tmpdir -t' ' -k2,2 -k3,3n -k4,4 -k5,5 \
	| awk '{print $2,$3,$3+1,$4,$5,$1}' \
	| tr " " "\t" > $tmpdir/$experiment.q30.umi.bed


###############################################################
# Filter out centromere, telomere, and blacklist regions
centbed=hg19-telomere-centromere-telomere.bed
blacklist=consensusBlacklist.bed
file=$tmpdir/$experiment.q30.umi.bed
bedtools intersect -v -a $file -b $centbed > $tmpdir/nocent.bed
bedtools intersect -v -a $tmpdir/nocent.bed -b $blacklist > $tmpdir/noblacklist.bed
mv $tmpdir/noblacklist.bed $tmpdir/$experiment.q30.umi.bed 

###############################################################
# Filter out PCR duplicates by getting unique UMIs
cat $tmpdir/$experiment.q30.umi.bed | cut -f-5 \
	| LC_ALL=C uniq -c \
	| awk '{print $2,$3,$4,$5,$6,$1}' \
	| tr " " "," > $tmpdir/$experiment.q30.umi.csv 
python umi_filtering.py $tmpdir/$experiment.q30.umi.csv $tmpdir/umi_filt1.txt
cat $tmpdir/umi_filt1.txt | cut -f1-3 | LC_ALL=C uniq -c \
	| awk '{OFS="\t"; print $2,$3,$4,$1}' > $tmpdir/q30_chr-loc-uniqueUMIcount.bed

#################################################################
# Get Summary for the data (number of reads after each filter, etc.) 
summaryfile=$tmpdir/summary.txt
echo "Number of fragments:" > $summaryfile
wc -l $tmpdir/r1oneline.fq >> $summaryfile
echo "Number of fragments with prefix:" >> $summaryfile
cat $tmpdir/*filtered* | paste - - | wc -l >> $summaryfile
echo "Alignment statistics:" >> $summaryfile 
samtools flagstat $tmpdir/$experiment.q30.sorted.bam >> $summaryfile
echo "Number of left and right cuts:" >> $summaryfile
cat $tmpdir/umi_filt1.txt | grep -v "_" | cut -f4 | sort | uniq -c >> $summaryfile
echo "Number of DSB locations:" >> $summaryfile
cat $tmpdir/q30_chr-loc-uniqueUMIcount.bed | grep -v "_" | wc -l >> $summaryfile


