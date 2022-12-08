# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 13:48:36 2022

@author: luisa
"""
import pandas as pd
import re

#count how many SNVs (rows are CDS) -> reused code snippet by Luis
def coding_genes (filename, outputfile):
    csv_file = pd.read_csv(filename)
    
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
    genes = {}
    SNVs_dict = {}
    for i in range(len(where3)):
        for j in range(len(where3[i])):
            if where3[i][j] == 'CDS' and SNVs_list2[i][j] == 'Nonsynonymous':
                genes[genename3[i][j]] = str(chrom[i])
                SNVs_dict[genename3[i][j]] = SNVs_list2[i][j]

    gene_name_values = [i for i in genes.keys()]
    chrom_values = [*genes.values()]
    SNVs_values = [*SNVs_dict.values()]

    lists = {'Gene': gene_name_values, 'Chrom': chrom_values, 'Change_Type': SNVs_values}
    df = pd.DataFrame(lists)
    df.to_csv(outputfile, index = False)
    
    return gene_name_values


CDS_normal_nsSNVs = coding_genes("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\Group3Normal.csv", "CDS_nsSNVs_normal.csv")
CDS_tumor_nsSNVs = coding_genes("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\Group3Tumor.csv", "CDS_nsSNVs_tumor.csv")


        