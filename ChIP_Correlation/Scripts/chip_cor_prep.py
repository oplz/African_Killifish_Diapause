#!/usr/env python3

k4 = {}
k27 = {}
pk27 = {}
rna = {}

with open("k4m3_avg.txt", 'r') as one, open("k27m3_avg.txt", 'r') as two, open("pk27m3.txt", 'r') as three, open("rna.txt", 'r') as four:
    for line in one:
        line = line.rstrip('\n').split('\t')
        if line[0] != "id":
            k4[line[0]] = line[1]
    
    for line in two:
        line = line.rstrip('\n').split('\t')
        if line[0] != "id":
            k27[line[0]] = line[1]
    
    for line in three:
        line = line.rstrip('\n').split('\t')
        if line[0] != "pK27me3_id":
            pk27[line[0]] = line[1]
    
    for line in four:
        line = line.rstrip('\n').split('\t')
        if line[0] != "RNA_id":
            rna[line[0]] = line[1]
            
with open("k4_vs_rna.txt", 'w') as out1, open("k27_vs_rna.txt", 'w') as out2, open("pk27_vs_rna.txt", 'w') as out3:
    out1.write("Gene\tChip_FC\tRNA_FC\n")
    out2.write("Gene\tChip_FC\tRNA_FC\n")
    out3.write("Gene\tChip_FC\tRNA_FC\n")
    
    for item in k4.keys():
        if item in rna.keys():
            out1.write(str(item) + '\t' + str(k4[item]) + '\t' + str(rna[item]) + '\n')
    
    for item in k27.keys():
        if item in rna.keys():
            out2.write(str(item) + '\t' + str(k27[item]) + '\t' + str(rna[item]) + '\n')
    
    for item in pk27.keys():
        if item in rna.keys():
            out3.write(str(item) + '\t' + str(pk27[item]) + '\t' + str(rna[item]) + '\n')
            
            