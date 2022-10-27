# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 23:04:30 2022

@author: luisa
"""
#SNP Nexus generator
import pandas as pd
import os
import csv
#1. Combine all OMI files into 1
#get all files
path = "C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\files"
file_list = os.listdir(path)
df_append = pd.DataFrame()
for file in file_list:
            df_temp = pd.read_csv(file)
            df_append = df_append.append(df_temp, ignore_index=True)
#export
df_append.to_csv('merged_omi.csv', index=False)

#2. Create a tab delimited text document containing the following
chrom_num = [] #store chrom num
pos_omi= [] #store position of variant
ref_base = [] #store ref base
alt_base = [] #store alternative base
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\files\\merged_omi.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            chrom_num.append(row[1])
            pos_omi.append(row[2])
            ref_base.append(row[5])
            alt_base.append(row[6])
        line_count+=1

num = [] # for each chrom store only the number
for i in range(len(chrom_num)):
    num.append(chrom_num[i][3])

#the word "chromosome" in the first col

#strand -> "1" for all since all alleles in the database are reported in the forward orientation
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\output_for_snpnexus.txt', 'w') as f:
    for i in range(len(num)):
        f.write("chromosome" + "\t" + num[i] + "\t" + pos_omi[i] + "\t" + ref_base[i] + "\t" + alt_base[i] + "\t" + "1" + "\n")

#Merge output with OMI file
dbsnp = []
pos_snpnexus = []
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\snp_nexus_outpu\\gen_coords_HW4.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            dbsnp.append(row[1])
            pos_snpnexus.append(row[3])
        line_count+=1
            
reordered_dbsnp = []
for i in range(len(pos_omi)):
    for j in range(len(pos_snpnexus)):
        if pos_snpnexus[j] == pos_omi[i]:
            reordered_dbsnp.append(dbsnp[j])
            break
        elif j == len(pos_snpnexus)-1 and pos_snpnexus[j] != pos_omi[i]:
            reordered_dbsnp.append("N/A")
        
pos_near_snpnexus = []            
overlapped_gene = []
annotation = []
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\snp_nexus_outpu\\near_gens_HW4.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            pos_near_snpnexus.append(row[2])
            overlapped_gene.append(row[3])
            annotation.append(row[5])
        line_count+=1

reordered_overlapped_gene = []
reordered_annotation = []
for i in range(len(pos_omi)):
    for j in range(len(pos_near_snpnexus)):
        if pos_near_snpnexus[j] == pos_omi[i]:
            reordered_overlapped_gene.append(overlapped_gene[j])
            reordered_annotation.append(annotation[j])
            break
        elif j == len(pos_near_snpnexus)-1 and pos_near_snpnexus[j] != pos_omi[i]:
            reordered_overlapped_gene.append("N/A")
            reordered_annotation.append("N/A")
            
#Add lists to csv to complete omi file



