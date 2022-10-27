# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 23:14:47 2022

@author: luisa
"""
import csv
import pandas as pd
#FATHMM-XF input generator
#1. Create a comma delimited text document containing the following (to use as input) #step 5
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

#using the lists, create the comma delimited file #step 6
with open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\output_for_fathmm.txt', 'w') as f:
    for i in range(len(num)):
        f.write(num[i] + "," + pos[i] + "," + ref_base[i] + "," + alt_base[i] + "\n")

#Merge output with OMI file #step 13
#create lists with needed information
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
        else: #create col headers
            coding_score.append("coding score")
            non_coding_score.append("non coding score")
            warning_msg.append("warning")
        line_count+=1       

#add cols to omi using pandas df #step 14
omi_file_df = pd.read_csv('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\merged_omi.csv') #CHANGE NAME TO FILE
omi_file_df["coding score"] = coding_score
omi_file_df['non coding score'] = non_coding_score
omi_file_df['warning'] = warning_msg

#convert df to csv #step 15
omi_file_df.to_csv('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 4\\pre_final_omi.csv', index=False) #CHANGE NAME TO FILE

"""
#add new cols other method DO NOT USE
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
"""