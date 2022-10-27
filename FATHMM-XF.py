# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 23:14:47 2022

@author: luisa
"""
import csv
#FATHMM-XF inout generator
#1. Create a comma delimited text document containing the following
chrom_num = [] #store chrom num
pos= [] #store position of variant
ref_base = [] #store ref base
alt_base = [] #store alternative base
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\files\\merged_omi.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            chrom_num.append(row[1])
            pos.append(row[2])
            ref_base.append(row[5])
            alt_base.append(row[6])
        line_count+=1

num = [] # for each chrom store only the number
for i in range(len(chrom_num)):
    num.append(chrom_num[i][3])

#the word "chromosome" in the first col

#strand -> "1" for all since all alleles in the database are reported in the forward orientation
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\output_for_fathmm.txt', 'w') as f:
    for i in range(len(num)):
        f.write(num[i] + "," + pos[i] + "," + ref_base[i] + "," + alt_base[i] + "\n")

#Merge output with OMI file
coding_score = []
non_coding_score = []
warning_msg = []
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\fathmm_output\\fathmm_hw4_output.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if line_count > 0:
            coding_score.append(row[4])
            non_coding_score.append(row[5])
            warning_msg.append(row[6])
        else:
            coding_score.append("coding score")
            non_coding_score.append("non coding score")
            warning_msg.append("warning")
        line_count+=1
        
        
#add new cols
from csv import writer
from csv import reader
def add_column_in_csv(input_file, output_file, transform_row):
    with open(input_file, 'r') as read_obj, \
                open(output_file, 'w', newline='') as write_obj:
            csv_reader = reader(read_obj)
            csv_writer = writer(write_obj)
            for row in csv_reader:
                transform_row(row, csv_reader.line_num)
                csv_writer.writerow(row)

add_column_in_csv('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi.csv', 'C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi1.csv', lambda row, line_num: row.append(coding_score[line_num - 1]))  
add_column_in_csv('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi1.csv', 'C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi2.csv', lambda row, line_num: row.append(non_coding_score[line_num - 1]))   
add_column_in_csv('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi2.csv', 'C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi3.csv', lambda row, line_num: row.append(warning_msg[line_num - 1]))         