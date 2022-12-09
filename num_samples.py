# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:45:06 2022

@author: lgraciamaz
"""

import re
import os

os.chdir("C:\\Users\\lgraciamaz\\Desktop\\Fall 2022\\Post-Genomic Analysis\\Project\\VCFFiles")
files = os.listdir()
vcf_files = []

for i in range(len(files)):
    if files[i].endswith(".vcf"):
        vcf_files += [files[i]]
        
regex = "##INDIVIDUAL=<NAME=TCGA-AB-(.+),"

Patient_ids = []
for filename in vcf_files:
    temp_file = open(filename, "r")
    file_content = temp_file.read()
    Patient_ids += [re.findall(regex, file_content)[0]]
    
def get_unique_numbers(numbers):
    list_of_unique_numbers = []
    unique_numbers = set(numbers)
    for number in unique_numbers:
        list_of_unique_numbers.append(number)
    return list_of_unique_numbers

unique_ids = get_unique_numbers(Patient_ids)

patient_occurrences = []
for i in range(len(unique_ids)):
    count = 0
    for j in range(len(Patient_ids)):
        if unique_ids[i] == Patient_ids[j]:
            count += 1
    patient_occurrences += [count]