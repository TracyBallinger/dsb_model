#!/bin/bash 

# make the tsvfile that has all of the SV for WGS cohorts. (Left out BOCA-UK and BRCA-UK because they aren't WGS.)
# gzip -dc structural/structural_somatic_mutation.BRCA-EU.tsv.gz | head -1 > tmph
# cat structural/structural_somatic_mutation.*.tsv.gz | gzip -dc | grep -v ^icgc_donor_id | cat tmph - | gzip > structural_somatic_mutation.WGS.tsv.gz

tsvfile=$1 # structural_somatic_mutation.WGS.tsv.gz
outdir=$2 # svseg_beds/panc_wgs
mkdir -p $outdir
svtype=any
gzip -dc $tsvfile | sed 1d | cut -f7,12,13,17,18,23 | \
	awk 'BEGIN{OFS="\t"; FS="\t"}{
	if ($2 != $4) print $2"\t"$3-1"\t"$3"\n"$4"\t"$5-1"\t"$5; 
	else if ($3>$5) print $2"\t"$5"\t"$3; 
	else print $2"\t"$3"\t"$5}' | \
	LC_ALL=C sort -k1,1 -k2,2n > $outdir/$svtype.bed

for svtype in deletion inversion insertion interchr intrachr translocation duplication; do 
	gzip -dc $tsvfile | sed 1d | cut -f7,12,13,17,18,23 | grep $svtype | \
	awk 'BEGIN{OFS="\t"; FS="\t"}{
	if ($2 != $4) print $2"\t"$3-1"\t"$3"\n"$4"\t"$5-1"\t"$5; 
	else if ($3>$5) print $2"\t"$5"\t"$3; 
	else print $2"\t"$3"\t"$5}' | \
	LC_ALL=C sort -k1,1 -k2,2n > $outdir/$svtype.bed
done

 
