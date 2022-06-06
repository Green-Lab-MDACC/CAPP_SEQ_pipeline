import sys
import pandas as pd
import os
import csv
import re 
import numpy as np
import math

sample_name = sys.argv[1]
VARIANTS_CALLED_TABLE = sys.argv[2]
FINAL_CALLS_TABLE = sys.argv[3]
PONPATH = sys.argv[4]
sample_id = sys.argv[5]
consensus = sys.argv[6]
PON = False
output_file_name = str(sample_name) + '.csv'
output_file_path = os.path.join(str(os.getcwd()), 'bnsugg', 'duplex_pipeline', 'data', sample_id, sample_name, consensus, output_file_name)
input_pon_file = PONPATH

samples = [sample_name]


j = 0

for i in range(0,len(samples)):

    input_file = FINAL_CALLS_TABLE
    df = pd.read_table(input_file)
    df.rename(columns={ df.columns[-1]: "DP" }, inplace = True)
    df.rename(columns={ df.columns[-2]: "AD" }, inplace = True)

    input_file = VARIANTS_CALLED_TABLE
    df2 = pd.read_table(input_file)
    df2.rename(columns={ df2.columns[-1]: "DP" }, inplace = True)
    df2.rename(columns={ df2.columns[-2]: "AD" }, inplace = True)

    df = pd.concat([df, df2])
    df.insert(0, 'SAMPLE', samples[i])
    
    #Delimit AD
    AD = df['AD'].str.split(',', expand=True)
    print(AD)
    df['RD'] = AD[0]
    df['VD'] = AD[1]
    #df['N'] = AD[2]
    df = df.drop_duplicates(subset=['POS'], keep='first')
    df = df.drop(['AD'],axis=1)

    if j > 0:
        og = pd.concat([og, df])
    else:
        og = df


    j = j + 1

df = og


if  PON is True:

    count = []
    PON = pd.read_excel(input_pon_file)
    pon_count = []
    PON = PON[(PON['DP'] > 100) & (PON['oVAF'] > 0.01)]

    for i in range(0,len(PON)):
        pos = PON.iloc[i]['POS']
        pon_count.append(len(PON[(PON['POS'] == pos)]))

    PON['Count'] = pon_count

    for i in range(0,len(df)):
        pos = df.iloc[i]['POS']
        g = PON[PON['POS'] == pos]['Count']
        if len(g) == 0:
            count.append(0)
        else:
            count.append(g.iloc[0])

    df["PoN Count"] = np.transpose(count)


df = df[(df['DP'] > 100)]
df.to_csv(output_file_path,index=False)




