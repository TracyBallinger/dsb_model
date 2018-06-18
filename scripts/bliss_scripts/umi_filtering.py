# -*- coding: utf-8 -*-

import sys
import csv
import numpy as np

filename = sys.argv[1]

with open(filename, 'rb') as f:
    reader = csv.reader(f)
    data = list(reader)

# sys.stderr.write("data is %d\n" % len(data))

# GROUP TOGETHER CONSECUTIVE IDENTICAL UMIS, ASSUMING THEIR PROXIMITY
# row_old = data[0]
# data_aggegated_by_umi_identity = [row_old]
# for row in data[1:]:
#    if (row_old[0]==row[0] and row_old[3]==row[3] and row_old[4]==row[4]):
#        row[5] = str(int(row[5])+int(row_old[5]))
#        if int(row_old[5]) > int(row_old[5]):
#            row[1:3] = row_old[1:3]
#        del data_aggegated_by_umi_identity[-1]
#    data_aggegated_by_umi_identity.append(row)
#    row_old = row

# GROUP TOGETHER CLOSE SPATIAL CONSECUTIVE READS WHOSE UMI DIFFERS AT MOST BY 2 MISMATCHES
data_aggegated_by_umi_similarity = []
skipped= []
space_gap = 30 
mm_gap = 2
rowi=1
oldi=0
while rowi < len(data):
	row=data[rowi]
	row_old=data[oldi]
	#sys.stderr.write("row_old: %s\n" % (row_old))
	#sys.stderr.write("row: %s\n" % (row))
	#sys.stderr.write("data_agg: %s\n" % str(data_aggegated_by_umi_similarity))
	s1 = row_old[4]
	s2 = row[4]
	numb_mismatches = sum(c1!=c2 for c1,c2 in zip(s1,s2))
	dist = abs(int(row[1])-int(row_old[1]))
	if (row_old[0]==row[0] and dist<=space_gap and row_old[3]==row[3] and numb_mismatches<=mm_gap):
		if int(row[5]) > int(row_old[5]):
			row_old[1:5] = row[1:5]
		row_old[5] = str(int(row[5])+int(row_old[5]))
		rowi=rowi+1
	else: 
		if dist > space_gap or row_old[0] != row[0]:
			data_aggegated_by_umi_similarity.append(row_old)
			if len(skipped) >0: 
				oldi=skipped.pop(0)
			else: 
				oldi=rowi
				rowi=rowi+1
		else: 
			skipped.append(rowi)
			rowi=rowi+1

row_old=data[oldi]
data_aggegated_by_umi_similarity.append(row_old)

thefile = open(sys.argv[2], 'wa') 
for item in data_aggegated_by_umi_similarity:
  thefile.write('\t'.join(item)+'\n')

print 'Done with filtering UMIs!'
