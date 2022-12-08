# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 00:22:51 2022

@author: luisa
"""
#DIFFERENCE: makes a df with all the information of each row for wach nsSNV
import pandas as pd
import re

#count how many SNVs (rows are CDS) -> reused code snippet by Luis
csv_file = pd.read_csv("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\Group3Tumor.csv")
new_df = pd.DataFrame(columns=csv_file.columns)

where = list(csv_file["Where"])
where2 = []
for i in range(len(where)):
    where2.append(re.sub('[^0-9a-zA-Z\s,]+', "", where[i]))
    
where3 = []
for i in range(len(where2)):
    where3.append(where2[i].split(","))
    
genename = list(csv_file['Gene_Name'])
genename2 = []
for i in range(len(genename)):
    genename2.append(re.sub('[^0-9a-zA-Z\s,]+', "", genename[i]))
genename3 = []
for i in range(len(genename2)):
    genename3.append(genename2[i].split(","))

change_type = list(csv_file['Change_Type'])
SNVs = []
for i in range(len(change_type)):
    SNVs.append(re.sub('[^0-9a-zA-Z\s,'']+', "", change_type[i])) #remove everything that is not letters
SNVs_list2 = []
for i in range(len(SNVs)):
    SNVs_list2.append(SNVs[i].split(","))
    
chrom = list(csv_file['chrom'])
unique_genes = []
genes = {}
SNVs_dict = {}
count = 0
for i in range(len(where3)):
    for j in range(len(where3[i])):
        if where3[i][j] == 'CDS' and SNVs_list2[i][j] == 'Nonsynonymous' and genename3[i][j] not in unique_genes:
            unique_genes.append(genename3[i][j])
            genes[genename3[i][j]] = str(chrom[i])
            SNVs_dict[genename3[i][j]] = SNVs_list2[i][j]
            new_df.loc[len(new_df.index)] = csv_file.iloc[i] #THIS IS DIFFERENT
            count +=1

gene_name_values = [i for i in genes.keys()]
chrom_values = [*genes.values()]
SNVs_values = [*SNVs_dict.values()]

lists = {'Gene': gene_name_values, 'Chrom': chrom_values, 'Change_Type': SNVs_values}
df = pd.DataFrame(lists)
df.to_csv("tumor_normal1.csv", index = False)

#new_df.to_csv("normal_nsSNVs_whole_row.csv", index = False)
new_df.to_csv("tumor_nsSNVs_whole_row.csv", index = False)
#return gene_name_values


#CDS_normal_nsSNVs = coding_genes("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\Group3Normal.csv", "test_normal1.csv")
#CDS_tumor_nsSNVs = coding_genes("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\Group3Tumor.csv", "test_tumor1.csv")
