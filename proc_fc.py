#! /usr/bin/env python

import sys

count_file = sys.argv[1]
summary_file = sys.argv[2]
gtf_file = sys.argv[3]
# output_file = sys.argv[3]

gene_id2gene_name = {}
with open(gtf_file, 'r') as hin:
    for line in hin:
        if line.startswith('#'): continue
        F = line.rstrip('\n').split('\t')
        
        if F[2] == "gene":
            infos = F[8].split(';')
            gene_id, gene_name = '', ''
            for elm in infos:
                elm = elm.strip(' ')
                if elm.startswith("gene_id"):
                    gene_id = elm.replace("gene_id ", '').strip('"')
                if elm.startswith("gene_name"):
                    gene_name = elm.replace("gene_name ", '').strip('"')
            gene_id2gene_name[gene_id] = gene_name
        

with open(summary_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] == "Assigned": 
            total_read_num = int(F[1])

total_read_ratio = float(total_read_num) / 1000000

gene2fpkm = {}
with open(count_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0].startswith('#'): continue
        if F[0] == "Geneid": continue
        gene_id = F[0]
        if gene_id not in gene_id2gene_name: continue
        gene = gene_id2gene_name[F[0]]
        gene_length = float(F[5])
        count = float(F[6])

        fpkm = count / ((gene_length / 1000) * total_read_ratio)
        if gene not in gene2fpkm: gene2fpkm[gene] = 0
        if fpkm > gene2fpkm[gene]: gene2fpkm[gene] = fpkm


for gene in gene2fpkm:
    print(gene + '\t' + str(round(gene2fpkm[gene], 4)))


