# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 19:11:39 2022

@author: lgraciamaz
"""
import re

normal_fasta = open("normal_genes_aa_sequences.fasta", 'r')
tumor_fasta = open("tumor_genes_aa_sequences.fasta", 'r')

normal_genes = normal_fasta.read()
tumor_genes = tumor_fasta.read()
normal_fasta.close()
tumor_fasta.close()

normal_sequence_only = re.findall(r'\n(.*)\n', normal_genes)
tumor_sequence_only = re.findall(r'\n(.*)\n', tumor_genes)

normal_amino_acid_occurrences = {"A" : 0,
                                 "G" : 0,
                                 "I" : 0,
                                 "L" : 0,
                                 "P" : 0,
                                 "V" : 0,
                                 "F" : 0,
                                 "W" : 0,
                                 "Y" : 0,
                                 "D" : 0,
                                 "E" : 0,
                                 "R" : 0,
                                 "H" : 0,
                                 "K" : 0,
                                 "S" : 0,
                                 "T" : 0,
                                 "C" : 0,
                                 "M" : 0,
                                 "N" : 0,
                                 "Q" : 0,
                                 "*" : 0}
normal_amino_acid_count = 0
for i in range(len(normal_sequence_only)):
    for j in range(len(normal_sequence_only[i])):
        normal_amino_acid_occurrences[normal_sequence_only[i][j]] += 1
        normal_amino_acid_count += 1
        
normal_amino_acid_frequencies = {"A" : normal_amino_acid_occurrences["A"] / normal_amino_acid_count,
                                 "G" : normal_amino_acid_occurrences["G"] / normal_amino_acid_count,
                                 "I" : normal_amino_acid_occurrences["I"] / normal_amino_acid_count,
                                 "L" : normal_amino_acid_occurrences["L"] / normal_amino_acid_count,
                                 "P" : normal_amino_acid_occurrences["P"] / normal_amino_acid_count,
                                 "V" : normal_amino_acid_occurrences["V"] / normal_amino_acid_count,
                                 "F" : normal_amino_acid_occurrences["F"] / normal_amino_acid_count,
                                 "W" : normal_amino_acid_occurrences["W"] / normal_amino_acid_count,
                                 "Y" : normal_amino_acid_occurrences["Y"] / normal_amino_acid_count,
                                 "D" : normal_amino_acid_occurrences["D"] / normal_amino_acid_count,
                                 "E" : normal_amino_acid_occurrences["E"] / normal_amino_acid_count,
                                 "R" : normal_amino_acid_occurrences["R"] / normal_amino_acid_count,
                                 "H" : normal_amino_acid_occurrences["H"] / normal_amino_acid_count,
                                 "K" : normal_amino_acid_occurrences["K"] / normal_amino_acid_count,
                                 "S" : normal_amino_acid_occurrences["S"] / normal_amino_acid_count,
                                 "T" : normal_amino_acid_occurrences["T"] / normal_amino_acid_count,
                                 "C" : normal_amino_acid_occurrences["C"] / normal_amino_acid_count,
                                 "M" : normal_amino_acid_occurrences["M"] / normal_amino_acid_count,
                                 "N" : normal_amino_acid_occurrences["N"] / normal_amino_acid_count,
                                 "Q" : normal_amino_acid_occurrences["Q"] / normal_amino_acid_count,
                                 "*" : normal_amino_acid_occurrences["*"] / normal_amino_acid_count}

tumor_amino_acid_occurrences = {"A" : 0,
                                 "G" : 0,
                                 "I" : 0,
                                 "L" : 0,
                                 "P" : 0,
                                 "V" : 0,
                                 "F" : 0,
                                 "W" : 0,
                                 "Y" : 0,
                                 "D" : 0,
                                 "E" : 0,
                                 "R" : 0,
                                 "H" : 0,
                                 "K" : 0,
                                 "S" : 0,
                                 "T" : 0,
                                 "C" : 0,
                                 "M" : 0,
                                 "N" : 0,
                                 "Q" : 0,
                                 "*" : 0}
tumor_amino_acid_count = 0
for i in range(len(tumor_sequence_only)):
    for j in range(len(tumor_sequence_only[i])):
        tumor_amino_acid_occurrences[tumor_sequence_only[i][j]] += 1
        tumor_amino_acid_count += 1
        
tumor_amino_acid_frequencies = {"A" : tumor_amino_acid_occurrences["A"] / tumor_amino_acid_count,
                                 "G" : tumor_amino_acid_occurrences["G"] / tumor_amino_acid_count,
                                 "I" : tumor_amino_acid_occurrences["I"] / tumor_amino_acid_count,
                                 "L" : tumor_amino_acid_occurrences["L"] / tumor_amino_acid_count,
                                 "P" : tumor_amino_acid_occurrences["P"] / tumor_amino_acid_count,
                                 "V" : tumor_amino_acid_occurrences["V"] / tumor_amino_acid_count,
                                 "F" : tumor_amino_acid_occurrences["F"] / tumor_amino_acid_count,
                                 "W" : tumor_amino_acid_occurrences["W"] / tumor_amino_acid_count,
                                 "Y" : tumor_amino_acid_occurrences["Y"] / tumor_amino_acid_count,
                                 "D" : tumor_amino_acid_occurrences["D"] / tumor_amino_acid_count,
                                 "E" : tumor_amino_acid_occurrences["E"] / tumor_amino_acid_count,
                                 "R" : tumor_amino_acid_occurrences["R"] / tumor_amino_acid_count,
                                 "H" : tumor_amino_acid_occurrences["H"] / tumor_amino_acid_count,
                                 "K" : tumor_amino_acid_occurrences["K"] / tumor_amino_acid_count,
                                 "S" : tumor_amino_acid_occurrences["S"] / tumor_amino_acid_count,
                                 "T" : tumor_amino_acid_occurrences["T"] / tumor_amino_acid_count,
                                 "C" : tumor_amino_acid_occurrences["C"] / tumor_amino_acid_count,
                                 "M" : tumor_amino_acid_occurrences["M"] / tumor_amino_acid_count,
                                 "N" : tumor_amino_acid_occurrences["N"] / tumor_amino_acid_count,
                                 "Q" : tumor_amino_acid_occurrences["Q"] / tumor_amino_acid_count,
                                 "*" : tumor_amino_acid_occurrences["*"] / tumor_amino_acid_count}