# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 01:42:39 2022

@author: luisa
"""

import csv
#FATHMM-XF input generator
#1. Create a comma delimited text document containing the following (to use as input) #step 5
chrom_num = [] #store chrom num
pos= [] #store position of variant
ref_base = [] #store ref base
alt_base = [] #store alternative base
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\normal_nsSNVs_whole_row.csv') as csv_file: #change for normal ans tumor
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            chrom_num.append(row[0])
            pos.append(row[1])
            ref_base.append(row[3])
            alt_base.append(row[5])
        line_count+=1

num = [] # for each chrom store only the number
for i in range(len(chrom_num)):
    if len(chrom_num[i]) > 4:
        num.append(chrom_num[i][-2:])
    else:
        num.append(chrom_num[i][3])

#using the lists, create the comma delimited file #step 6
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\inputs\\input_for_fathm_normal.txt', 'w') as f:
    for i in range(len(num)):
        f.write(num[i] + "," + pos[i] + "," + ref_base[i] + "," + alt_base[i] + "\n")
