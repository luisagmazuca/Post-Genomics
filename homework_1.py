# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:56:06 2022

@author: luisa
"""
import vcf

vcf_reader = vcf.Reader(open('C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 1\\example_VCF.vcf', 'r'))
records = []
for record in vcf_reader:
    print(record) #Print all records
    records.append(record) 
   
#Information for the first record in the vcf
print("CHROM:", records[0].CHROM)
print("POS:", records[0].POS)
print("ID:", records[0].ID)
print("REF", records[0].REF)
print("ALT", records[0].ALT)
print("QUAL:", records[0].QUAL)
print("FILTER:", records[0].FILTER)
print("INFO:", records[0].INFO)
print("samples:", records[0].samples)
print("genotype:", records[0].genotype)

#Meta information regarding the vcf
print("file date (yyymmdd):", vcf_reader.metadata['fileDate'])
print("center:", vcf_reader.metadata['center'])
print("reference:", vcf_reader.metadata['reference'])
#vcf_reader.metadata #for all metadata
print(vcf_reader.filters)
print(vcf_reader.formats)
print(vcf_reader.infos)

#number of samples
print("There are", len(vcf_reader.samples), "samples")
#sample names
print("The samples are:", vcf_reader.samples)
#number of heterozygous genotypes
print("There are", len(records), "heterozygous genotypes")
#variant type
print("The variant type is", record.var_type)






