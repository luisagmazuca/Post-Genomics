# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 18:11:01 2022

@author: luisa
"""
import csv

#2. Create a tab delimited text document containing the following (to use as input) #step 3
chrom_num = [] #store chrom num
pos_omi= [] #store position of variant
ref_base = [] #store ref base
alt_base = [] #store alternative base
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\tumor_nsSNVs_whole_row.csv') as csv_file: #change file name
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            chrom_num.append(row[0]) #'chrom'
            pos_omi.append(row[1]) #'left'
            ref_base.append(row[5]) #'var_seq1'
            alt_base.append(row[6]) #'var_seq2'
        line_count+=1

num = [] # for each chrom store only the number of the chrom
for i in range(len(chrom_num)):
    if len(chrom_num[i]) > 4:
        num.append(chrom_num[i][-2:])
    else:
        num.append(chrom_num[i][3])

#using the lists, create the tab delimited file #step 4
#put the word "chromosome" in the first col 
#in the last col strand -> "1" for all since all alleles in the database are reported in the forward orientation
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\inputs\\input_for_snpnexus_tumor.txt', 'w') as f: #change name of file for tumor and normal
    for i in range(len(num)):
        f.write("chromosome" + "\t" + num[i] + "\t" + pos_omi[i] + "\t" + ref_base[i] + "\t" + alt_base[i] + "\t" + "1" + "\n")
