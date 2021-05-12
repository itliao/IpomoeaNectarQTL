#!/usr/bin/env python
# coding: utf-8

#python conversion from interactive use in jupyter notebook

# Goal - extract variance values to calculate broad-sense heritability
# and to put into tabular form the mean, standard deviation, and standard error
# for C (cordatotriloba), L (lacunosa), and F5 RILs

import pandas as pd
import numpy as np
import random
import re
import os
import math


# working directory
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/F5_RILs_data/DataAnalysis/CleanedFiles/3R_final3_200803")
path = os.getcwd()
path


# list files in working directory
ls

############################
# calculating heritability #
############################

# create a list to put the variance values
varList = []


# open file and loop through file output from R and input values in a list
f = open("ANOVA_hert_variance_3R_final.txt")
for line in f:
    words = line.split(" ")
    #print words
    if words[0] == '[1]':
        key = words[1]
        print(key)
        varList.append(key.split("\n")[0])
    elif words[0] == "":
        if words[1] == 'F5line':
            var1 = words[5]
            print(var1)
            varList.append(var1)
        elif words[1] == 'Residual':
            var2 = words[14]
            print(var2)
            varList.append(var2)
varList      


# reformat list into list of lists
newList = []

i=0
while i < len(varList):
    newList.append([varList[i], varList[i+1], varList[i+2]])
    #print newList
    i = i+3
newList


#create dataframe from list
varDF = pd.DataFrame.from_records(newList)


#rename columns
varDF.columns=["Trait","F5lineVariance","ResVariance"]
varDF


#change values in columns to floats
varDF["F5lineVariance"] = varDF["F5lineVariance"].astype(float)
varDF["ResVariance"] = varDF["ResVariance"].astype(float)


# calculate the total variance and the broad-sense heritability
varDF["TotalVar"] = varDF["F5lineVariance"] + varDF["ResVariance"]
varDF["Heritability"] = varDF["F5lineVariance"]/varDF["TotalVar"]
varDF


############################################
# mean, standard deviation, standard error #
############################################

# create new list and dictionary of summary statistics
sumList = []
sumDict = {}


# open and loop through file of summary statistics
f = open("F5_3R_summary_final.txt")
for line in f:
    words = line.split(" ")
    #print words
    if words[0] == '[1]':
        key = words[1].split("\n")
        key2 = key[0]
    elif words[0] in ["1","2","3"]:
        print(words)
        value = words[-1].split("\n")
        print(value[0])
        if key2 in sumDict:
            sumDict[key2].append(value[0])
        else:
            sumDict[key2]=[value[0]]
sumDict


# convert dictionary to dataframe
sumDF = pd.DataFrame.from_dict(sumDict)


#transpose dataframe
sumDFt = sumDF.transpose()
sumDFt


#rename columns
sumDFt = sumDFt.reset_index()
#rename columns
sumDFt = sumDFt.rename(columns={0:'C_mean', 1:'F5_mean', 2:'L_mean', 3:'C_sd', 4:'F5_sd', 5:'L_sd', 6:"C_se", 7:'F5_se', 8:"L_se", 9:'C_se2', 10:'F5_se2', 11:'L_se2'})
sumDFt


##################################################
# merge heritability and summary statistics file #
##################################################

mergeHsum = sumDFt.merge(varDF, left_on="index", right_on="Trait", how="outer")
mergeHsum


#write to file
mergeHsum.to_csv("SumStats_Heritability_3R_200803.csv", sep = ",", encoding='utf-8', index=False)


