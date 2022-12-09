# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:54:47 2022

@author: lgraciamaz
"""

import pandas as pd
import re

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
        
    chrom = list(csv_file['chrom'])
    genes = {}
    for i in range(len(where3)):
        for j in range(len(where3[i])):
            if where3[i][j] == 'CDS':
                genes[genename3[i][j]] = str(chrom[i])
    
    df = pd.DataFrame(genes.items(), columns=['Gene', 'Chrom'])
    df.to_csv(outputfile, index = False)


#gene_names = coding_genes("C:\\Users\\lgraciamaz\\Desktop\\Fall_2022\\Post-Genomic Analysis\\Project\\VCFFiles\\output_csvs\\fe0333d9-1c86-437d-813f-980313903f35_tumor.csv")

gene_names_big_file_normal = coding_genes("C:\\Users\\lgraciamaz\\Desktop\\Fall_2022\\Post-Genomic Analysis\\Project\\Group3_normal.csv", "normal_genes.csv")

gene_names_big_file_tumor = coding_genes("C:\\Users\\lgraciamaz\\Desktop\\Fall_2022\\Post-Genomic Analysis\\Project\\Group3_tumor.csv", "tumor_genes.csv")


