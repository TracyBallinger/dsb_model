#!/bin/bash 


echo -e "cohort\tTot\tWGS\tWGA\tWXS\tdeletion\tinversion\tinsertion\ttand_dup\tinterchr\tintrachr\ttranslocation" > counts_per_cohort_sv.txt 


for tsvfile in `ls structural/structural_somatic_mutation.*.tsv.gz`; do 
	c=`echo $tsvfile | awk '{split($1, a, "."); print a[2]}'`
	nsamp=`gzip -dc $tsvfile | sed 1d | cut -f4 | sort -u | wc -l`
	wga=`gzip -dc $tsvfile | sed 1d | cut -f4,23 | grep WGA | cut -f1 | sort -u | wc -l`
	wgs=`gzip -dc $tsvfile | sed 1d | cut -f4,23 | grep WGS | cut -f1 | sort -u | wc -l`
	wxs=`gzip -dc $tsvfile | sed 1d | cut -f4,23 | grep WXS | cut -f1 | sort -u | wc -l`
	# get the types of mutations
	gzip -dc $tsvfile | cut -f7 | sort | uniq -c > tmpmutype
	del=`grep "deletion" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $del ]] && del=0
	inv=`grep "inversion" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $inv ]] && inv=0
	ins=`grep "insertion" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $ins ]] && ins=0
	dup=`grep "tandem duplication" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $dup ]] && dup=0
	inter=`grep "interchromosomal" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $inter ]] && inter=0
	intra=`grep "intrachromosomal" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $intra ]] && intra=0
	tran=`grep "translocation" tmpmutype | awk '{x+=$1}END{print x}'`
	[[ -z $tran ]] && tran=0
	echo -e "$c\t$nsamp\t$wgs\t$wga\t$wxs\t$del\t$inv\t$ins\t$dup\t$inter\t$intra\t$tran" >> counts_per_cohort_sv.txt

done 
