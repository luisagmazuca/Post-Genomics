# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:45:17 2022

@author: luisa
"""

#SECTION 6C
import pandas as pd
from numpy import log as ln

#get lenghts for each gene
csv_file = pd.read_csv("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\prelim_file2_tumor.csv")
gene = list(csv_file['Gene name'])

#get lenghths from file
csv_file_seq_len = pd.read_csv("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\CDS_AML_Tumor_Sequences_Lengths.csv")
all_genes = list(csv_file_seq_len['Gene'])
gene_len = list(csv_file_seq_len['Length'])

#return indexes after searching
len_list = []
for i in range(len(gene)):
    if gene[i] not in all_genes:
        len_list.append('Not found')
    else:
        temp_index = all_genes.index(gene[i])
        len_list.append(gene_len[temp_index])
    
#append new col
csv_file['Seq. length'] = len_list

#apply formula
#get counts from opposite variant type
csv_file_normal = pd.read_csv("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\prelim_file2_normal.csv") #change variable name

#get pathogenic variants
#csv_file_pathogenic = pd.read_csv("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Project\\pathogenic_nsSNVs_tumor.csv")

#function for formula
def scoring_formula(i):
    curr_v = list(csv_file['Position'])[i]
    #score = 0
    find_in_normal = list(csv_file_normal['Position']) #change variable name
    if curr_v not in find_in_normal: #change variable name
        s_normal = 0
    else:
        temp_index = find_in_normal.index(curr_v) #change variable name
        s_normal = float(list(csv_file_normal['No. Patients'])[temp_index]) #change variable name
    s_tumor = list(csv_file['No. Patients'])[i] #change variable name
    if list(csv_file['FATHMM Score'])[i] == '--' or list(csv_file['Seq. length'])[i] == 'Not found':
        score = 'not calculated'
    else:
        f_score = float(list(csv_file['FATHMM Score'])[i]) #some do not have fathmm score
        s_len = float(ln(list(csv_file['Seq. length'])[i])) #some do not have length
        score = (1/s_len)*(f_score *(s_tumor - s_normal))
    return score
    
#call formula in loop to return list and append to file
cal_score = []
for i in range(len(gene)):
    curr_score = scoring_formula(i)
    cal_score.append(curr_score)
    
#append col to df
csv_file['Scoring func.'] = cal_score

#get rows with calculated score only
new_df = pd.DataFrame(columns=csv_file.columns)
for i in range(len(gene)):
    if list(csv_file['Scoring func.'])[i] != 'not calculated':
        new_df.loc[len(new_df.index)] = csv_file.iloc[i]

#rearrange highest to lowest
test = new_df.sort_values(['Scoring func.'], ascending=[False])
test.to_csv("top_tumor.csv", index = False)

#Export all cols as a csv file        
csv_file.to_csv("scoring_func_tumor.csv", index = False)



