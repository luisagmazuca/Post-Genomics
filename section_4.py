# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:41:43 2022

@author: lgraciamaz
"""
import pandas as pd
import re

normal = pd.read_csv("Group3_normal.csv")[["ref_seq","alt_seq", "Where"]].copy()
tumor = pd.read_csv("Group3_tumor.csv")[["ref_seq","alt_seq", "Where"]].copy()

where = list(normal["Where"])
where2 = []
for i in range(len(where)):
    where2.append(re.sub('[^0-9a-zA-Z\s,]+', "", where[i]))
where3 = []
for i in range(len(where2)):
    where3.append(where2[i].split(", "))

indexes_to_delete_normal = []
for i in range(len(where3)):
    count = 0
    for j in range(len(where3[i])):
        if where3[i][j] == 'CDS':
            count += 1
    if count == 0:
        indexes_to_delete_normal.append(i)
            
normal = normal.drop(labels = indexes_to_delete_normal, axis = 0)
normal = normal.reset_index()


twhere = list(tumor["Where"])
twhere2 = []
for i in range(len(twhere)):
    twhere2.append(re.sub('[^0-9a-zA-Z\s,]+', "", twhere[i]))
twhere3 = []
for i in range(len(twhere2)):
    twhere3.append(twhere2[i].split(", "))
    
for i in range(len(tumor)):
    tumor["Where"][i] = twhere3[i]

indexes_to_delete_tumor = []
for i in range(len(twhere3)):
    count = 0
    for j in range(len(twhere3[i])):
        if twhere3[i][j] == 'CDS':
            count += 1
    if count == 0:
        indexes_to_delete_tumor.append(i)
            
tumor = tumor.drop(labels = indexes_to_delete_tumor, axis = 0)
tumor = tumor.reset_index()
occurrences_normal = {"A->T" : 0,
                      "A->C" : 0,
                      "A->G" : 0,
                      "T->A" : 0,
                      "T->C" : 0,
                      "T->G" : 0,
                      "C->A" : 0,
                      "C->T" : 0,
                      "C->G" : 0,
                      "G->A" : 0,
                      "G->T" : 0,
                      "G->C" : 0}

frequencies_normal = {"A->T" : 0,
                      "A->C" : 0,
                      "A->G" : 0,
                      "T->A" : 0,
                      "T->C" : 0,
                      "T->G" : 0,
                      "C->A" : 0,
                      "C->T" : 0,
                      "C->G" : 0,
                      "G->A" : 0,
                      "G->T" : 0,
                      "G->C" : 0}

observed_bases_normal = {"A" : 0,
                         "T" : 0,
                         "C" : 0,
                         "G" : 0}

for i in range(len(normal)):
    occurrences_normal[str(normal["ref_seq"][i])+"->"+str(normal["alt_seq"][i])] += 1
    observed_bases_normal[str(normal["ref_seq"][i])] += 1
    
for key in frequencies_normal:
    frequencies_normal[key] = occurrences_normal[key] / len(normal)
    
occurrences_tumor = {"A->T" : 0,
                      "A->C" : 0,
                      "A->G" : 0,
                      "T->A" : 0,
                      "T->C" : 0,
                      "T->G" : 0,
                      "C->A" : 0,
                      "C->T" : 0,
                      "C->G" : 0,
                      "G->A" : 0,
                      "G->T" : 0,
                      "G->C" : 0}

frequencies_tumor = {"A->T" : 0,
                      "A->C" : 0,
                      "A->G" : 0,
                      "T->A" : 0,
                      "T->C" : 0,
                      "T->G" : 0,
                      "C->A" : 0,
                      "C->T" : 0,
                      "C->G" : 0,
                      "G->A" : 0,
                      "G->T" : 0,
                      "G->C" : 0}

observed_bases_tumor = {"A" : 0,
                         "T" : 0,
                         "C" : 0,
                         "G" : 0}


for i in range(len(tumor)):
    occurrences_tumor[str(tumor["ref_seq"][i])+"->"+str(tumor["alt_seq"][i])] += 1
    observed_bases_tumor[str(tumor["ref_seq"][i])] += 1
    
for key in frequencies_tumor:
    frequencies_tumor[key] = occurrences_tumor[key] / len(tumor)
    
mutation_levels_normal = {"A->T" : occurrences_normal["A->T"] / observed_bases_normal["A"],
                      "A->C" : occurrences_normal["A->C"] / observed_bases_normal["A"],
                      "A->G" : occurrences_normal["A->G"] / observed_bases_normal["A"],
                      "T->A" : occurrences_normal["T->A"] / observed_bases_normal["T"],
                      "T->C" : occurrences_normal["T->C"] / observed_bases_normal["T"],
                      "T->G" : occurrences_normal["T->G"] / observed_bases_normal["T"],
                      "C->A" : occurrences_normal["C->A"] / observed_bases_normal["C"],
                      "C->T" : occurrences_normal["C->T"] / observed_bases_normal["C"],
                      "C->G" : occurrences_normal["C->G"] / observed_bases_normal["C"],
                      "G->A" : occurrences_normal["G->A"] / observed_bases_normal["G"],
                      "G->T" : occurrences_normal["G->T"] / observed_bases_normal["G"],
                      "G->C" : occurrences_normal["G->C"] / observed_bases_normal["G"]}

mutation_levels_tumor = {"A->T" : occurrences_tumor["A->T"] / observed_bases_tumor["A"],
                      "A->C" : occurrences_tumor["A->C"] / observed_bases_tumor["A"],
                      "A->G" : occurrences_tumor["A->G"] / observed_bases_tumor["A"],
                      "T->A" : occurrences_tumor["T->A"] / observed_bases_tumor["T"],
                      "T->C" : occurrences_tumor["T->C"] / observed_bases_tumor["T"],
                      "T->G" : occurrences_tumor["T->G"] / observed_bases_tumor["T"],
                      "C->A" : occurrences_tumor["C->A"] / observed_bases_tumor["C"],
                      "C->T" : occurrences_tumor["C->T"] / observed_bases_tumor["C"],
                      "C->G" : occurrences_tumor["C->G"] / observed_bases_tumor["C"],
                      "G->A" : occurrences_tumor["G->A"] / observed_bases_tumor["G"],
                      "G->T" : occurrences_tumor["G->T"] / observed_bases_tumor["G"],
                      "G->C" : occurrences_tumor["G->C"] / observed_bases_tumor["G"]}