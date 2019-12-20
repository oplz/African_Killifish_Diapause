#!/usr/env python3

import sys

input = sys.argv[1]
output = sys.argv[2]
refer = sys.argv[3]
gene_2_fold = {}
sig_peaks = []

with open(input, 'r') as inner, open(output, 'w') as outer, open(refer, 'r') as ref:
	for line in ref:
		line = line.rstrip('\n').split('\t')
		if line[0] != 'id':
			sig_peaks.append(line[0])
	
	for line in inner:
		line = line.rstrip('\n').split('\t')
		if line[0] != 'seqnames':
			gene = line[17]
			fold = line[8]
			
			if gene in gene_2_fold.keys():
				gene_2_fold[gene].append(fold)
			else:
				gene_2_fold[gene] = []
				gene_2_fold[gene].append(fold)
	
	outer.write('id\tfc\n')
	for entry in gene_2_fold.keys():
		aver = 0
		leng = len(gene_2_fold[entry])
		for item in gene_2_fold[entry]:
			if aver == 0:
				aver = float(item)
			else: 
				aver = aver + float(item)
		aver = str(aver/leng)
		if entry in sig_peaks:
			outer.write(entry + '\t' + aver + '\n')
		
		
			