#!/usr/bin/env python
# coding: utf-8

#python conversion from interactive use in jupyter notebook

#Goal - to randomly choose 3 measurements in individuals with 3+ phenotypic measurements
#to prevent bias in collection methods

#randomizing script: https://www.geeksforgeeks.org/how-to-randomly-select-rows-from-pandas-dataframe/

import pandas as pd
import numpy as np
import random
import re
import os

#path to working directory
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/F5_RILs_data/DataAnalysis/CleanedFiles/")

#read in files using pandas
#df_f = flower data, only individuals with at least 3 flowers phenotyped
#df_n = nectar data, only individuals with at least 3 flowers phenotyped
df_f = pd.read_csv('F5Flower_Mod1-3_rm3less_190905.txt', sep='\t', engine='python', index_col=None)
df_n = pd.read_csv('F5Nectar_rm3less_190905.txt', sep='\t', engine='python', index_col=None)

#check if file imported correctly
df_f
df_n

#for df_f = individual in column "3" (start counting from 0)
#for df_n = individual in column "4"
print(df_f.loc[[455]])
print(df_f.iloc[455,3])
print(df_n.loc[[455]])
print(df_n.iloc[455,4])

#create dictionary of lists with the "key" as the individual plant
#len(df2.index) == number of rows
#df2.iloc[a,3] == column 3 is the Ind, a refers to the row
#append each dictionary list with the entire row

#for df_f == individual in column "3" (start counting from 0)
DictF = {}
for a in range(len(df_f.index)):
    key = df_f.iloc[a,3]
    if key in DictF:
        DictF[key].append(df_f.loc[[a]])
    else:
        DictF[key] = [df_f.loc[[a]]]
        
#for df_n == individual in column "4"        
DictN = {}
for a in range(len(df_n.index)):
    key = df_n.iloc[a,4]
    if key in DictN:
        DictN[key].append(df_n.loc[[a]])
    else:
        DictN[key] = [df_n.loc[[a]]] 

# check each dictionary - if call the individual number, what the output is
print(DictF['163'])
print(DictN['163'])

#can only sample 3; if less than 3 samples, cannot randomly sample 3
random.sample(DictF["163"], k=3)

#df1.sample(n = 3, replace = False) - use only for data frames
#because created a dictionary of LISTS, need to use random for lists, random.sample()
#Used for random sampling without replacement. (https://www.geeksforgeeks.org/python-random-sample-function/)
#https://pynative.com/python-random-choice/

#create new dataframes append new data
dfF3=pd.DataFrame()
dfN3=pd.DataFrame()

for key in DictF:
    dfF3.append(random.sample(DictF[str(key)], k=3))
    dfF3=dfF3.append(random.sample(DictF[str(key)], k=3))

for key in DictN:
    dfN3.append(random.sample(DictN[str(key)], k=3))
    dfN3=dfN3.append(random.sample(DictN[str(key)], k=3))
    
#check the output of the newly made dataframes
dfF3
dfN3

#write dataframes to tab-delimited file
dfF3.to_csv('F5random3_Flower_190905.txt', sep='\t', encoding='utf-8', index=False)
dfN3.to_csv('F5random3_Nectar_190905.txt', sep='\t', encoding='utf-8', index=False)


