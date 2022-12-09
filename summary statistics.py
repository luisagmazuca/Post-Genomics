import pandas as pd
normal = pd.read_csv("/Users/Michael Orosz/OneDrive/Desktop/Final Project/normal.csv")
tumor = pd.read_csv("/Users/Michael Orosz/OneDrive/Desktop/Final Project/tumor.csv")
PatientID = list(normal["Patient_ID"])
PatientID1 = list(tumor["Patient_ID"])
uniqueID = set(PatientID)
uniqueID1 = set(PatientID1)

for element in uniqueID:
    normal.insert(normal.shape[1], str(element), 0)
for element1 in uniqueID1:
    tumor.insert(tumor.shape[1], str(element1), 0)
for i in range(normal.shape[0]):
    temp = normal["Patient_ID"][i]
    normal[str(temp)][i] = 1
for j in range(tumor.shape[0]):
    temp1 = tumor["Patient_ID"][j]
    tumor[str(temp1)][j] = 1
grouped_normal = normal.groupby(['chrom','left','ref_seq','alt_seq','var_seq1','var_seq2','Gene_Name']).agg('sum').reset_index()
grouped_tumor = tumor.groupby(['chrom','left','ref_seq','alt_seq','var_seq1','var_seq2','Gene_Name']).agg('sum').reset_index()
grouped_normal = grouped_normal.drop(["Unnamed: 0", "Unnamed: 0.1", "count1","count2"], axis=1)
normal=grouped_normal.iloc[:, 9:64]
normT = normal.T
normplot = normT[['Normal SNV']]
tumor=grouped_tumor.iloc[:, 9:64]
tumor.loc['Tumor SNV',:]= tumor.sum(axis=0)
tumorT = tumor.T
tumorplot = tumorT[['Tumor SNV']]
joined_df=tumorplot.merge(normplot, how='left', left_index=True,right_index=True)
import numpy as np
joined_df['Log(Tumor SNV)'] = np.log(joined_df['Tumor SNV'])
joined_df['Log(Normal SNV)'] = np.log(joined_df['Normal SNV'])
logplot = joined_df[['Log(Tumor SNV)','Log(Normal SNV)']]
logsortplot = logplot.sort_values(["Log(Tumor SNV)","Log(Normal SNV)"], axis = 0, ascending = True)
import matplotlib.pyplot as plt

plt.plot(logsortplot,'.')
plt.title('Tumor vs Normal')
plt.ylabel('Log(SNV Count)')
plt.tick_params(
    axis='x',          # ticks on x axis are lame
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
location = 0 # For the best location
legend_drawn_flag = True
plt.legend(["Tumor", "Normal"], loc=0, frameon=legend_drawn_flag)