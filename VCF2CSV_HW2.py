import os

import pandas as pd
# To import vcf you may need to manually add "PyVCF" package in the preferences settings > Project: > Python Interpreter
import vcf


def writecsv(inputfile):
    vcf_reader = vcf.Reader(open(inputfile,'r'))  # read vcf file
    #use nested list to store variant information, one for tumor, one for normal
    normal_var=[]
    tumor_var=[]
    var_score = [28]
    refct = []
    altct = []
    refct1 = []
    altct1 = []

    for record in vcf_reader:
        curr_var = []
        left = []
        right = []
        ref_seq = []
        alt_seq = []

        tumorAD = record.samples[1]["AD"]
        normalAD = record.samples[0]["AD"]
        ########################################### tumor_var #############################################
        #######check for 2 ALT options:
        if len(tumorAD) == 3:
            sum = tumorAD[0] + tumorAD[1] + tumorAD[2]
            x = tumorAD[0]
            y = tumorAD[1]
            z = tumorAD[2]
            # if two out of the three are not greater than 0.05 throw away variant

            ## insert conditional statement that uses the defined variables x, y, z, and sum
            # that passes over variant if the above is not true.
            if (x < 0.05 and y < 0.05 or x < 0.05 and z < 0.05 or y < 0.05 and z < 0.05):
                pass






        ########only 1 base for the ALT sequence
        else:
            sum = tumorAD[0] + tumorAD[1]
            refcount = tumorAD[0]
            altcount = tumorAD[1]

            # both the REF and ALT have to be greater
            if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
                curr_var.append(record.CHROM)
                left.append(record.POS)
                right.append(record.POS + len(record.ALT))
                ref_seq.append(record.REF)
                alt_seq.append(str((record.ALT[0])))
                refct.append(refcount)
                altct.append(altcount)

                lists = curr_var + left + right + ref_seq + alt_seq + var_score
                tumor_var.append(lists)


            # throw away variant if not
            else:
                pass

        ############################# normal_var ##############################3
        curr_var1 = []
        left1 = []
        right1 = []
        ref_seq1 = []
        alt_seq1 = []
        if len(normalAD) == 3:
            sum = normalAD[0] + normalAD[1] + normalAD[2]
            x = normalAD[0]
            y = normalAD[1]
            z = normalAD[2]
            # if two out of the three are not greater than 0.05 throw away variant

            ## insert conditional statement that uses the defined variables x, y, z, and sum
            # that passes over variant if the above is not true.
            if (x < 0.05 and y < 0.05 or x < 0.05 and z < 0.05 or y < 0.05 and z < 0.05):
                pass





        ########only 1 base for the ALT sequence
        else:
            sum = normalAD[0] + normalAD[1]
            refcount = normalAD[0]
            altcount = normalAD[1]

            # both the REF and ALT have to be greater
            if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
                curr_var1.append(record.CHROM)
                left1.append(record.POS)
                right1.append(record.POS + len(record.ALT))
                ref_seq1.append(record.REF)
                alt_seq1.append(str((record.ALT[0])))
                refct1.append(refcount)
                altct1.append(altcount)

                lists = curr_var1 + left1 + right1 + ref_seq1 + alt_seq1 + var_score
                normal_var.append(lists)


            # throw away variant if not
            else:
                pass


    #transfer list into dat zaframe, it would be easier for following manipulation
    normal=pd.DataFrame(normal_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score'])
    tumor=pd.DataFrame(tumor_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score'])

    normal_seq2 = []
    tumor_seq2 = []

    normal.insert(5, "count1", altct1)
    normal.insert(6, "count2", refct1)
    for i in range(normal.shape[0]):

        if normal['ref_seq'][i] > normal['var_seq1'][i]:
            normal_seq2.append(normal['ref_seq'][i])

        else:
            normal_seq2.append(normal['var_seq1'][i])
            normal['var_seq1'][i] = normal['ref_seq'][i]

            temp = normal['count1'][i]
            normal['count1'][i] = normal['count2'][i]
            normal['count2'][i] = temp

    normal.insert(5, "var_seq2", normal_seq2)

    tumor.insert(5, "count1", altct)
    tumor.insert(6, "count2", refct)
    for j in range(tumor.shape[0]):
        if tumor['ref_seq'][j] > tumor['var_seq1'][j]:
            tumor_seq2.append(tumor['ref_seq'][j])

        else:
            tumor_seq2.append(tumor['var_seq1'][j])
            tumor['var_seq1'][j] = tumor['ref_seq'][j]
            # switching count values because the ref order is now var_seq1
            temp = tumor['count1'][j]
            tumor['count1'][j] = tumor['count2'][j]
            tumor['count2'][j] = temp

    tumor.insert(5, "var_seq2", tumor_seq2)

    #inserts and sets the var_index using a list of number within the range of the rows
    tumor.insert(0, "var_index", list(range(tumor.shape[0])))
    normal.insert(0, "var_index", list(range(normal.shape[0])))



    fil = inputfile.split(".")
    id = fil[0]

    tumor.to_csv(id + "_tumor.csv", index=False)
    normal.to_csv(id + "_normal.csv", index=False)
    return



if __name__=='__main__':
    # os.chdir("/Users/abataycan/Documents/PhD/Fall 2022/Post Genomics/Lab 1/VCF Files")
    # files = os.listdir()
    # for file in files:
    #     writecsv(file)

    ## Insert a code that reads/calls the given vcf files into the defined functioned
    # above is an example of how this can be done.
    os.chdir("C:\\Users\\luisa\\Documents\\Fall 2022\\Post Gen Analysis\\Homework 2\\vcf_files")
    files = os.listdir()
    for file in files:
         writecsv(file)
    print("Done")
