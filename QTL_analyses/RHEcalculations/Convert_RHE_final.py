#!/usr/bin/env python
# coding: utf-8

#python conversion from interactive use in jupyter notebook

# Goals: 
# (1) to convert text output from R of the number of each genotype, the mean phenotypic
# value into a dataframe, then calculate the RHE (relative homozygous effect)
# (2) merge information with what the 11/22 homozygous genotypes actually are associated 
# with (CC/LL) to determine actual direction of RHE

import pandas as pd
import random
import re
import os
import math


# working directory
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles")


# path to working directory
path = os.getcwd()
path


# read in file of parental differences of each phenotype
parent = pd.read_csv(path+"/200817_RQTL2/Parent_Means_200821.txt", sep="\t",engine='python') 
parent


###############################
# GWS (conservative set) QTLs #
###############################


#extract information from text file of the n, mean phenotypic value at peak, etc.
#before running code, double check text files 
#needed to modify to depending on the number outputs and remove the space (if 2 digits instead of 3 digits)
#and whether 22 was placed first

meanList = []
f = open("200805_RQTL/jitter/Final_GenoPhenoMeansNum_qtlCon_200827m.txt")
for line in f:
    words = line.split(" ")
    #print words
    if words[0] == '[1]':
        if words[1][0] == '"':
            key = words[1].split("\n")
            key2 = key[0]
            #print key2
            meanList.append(key2)
        else:
            value = words[1].split("\n")
            #print value[0]
            meanList.append(value[0])
    elif words[0] == "11":
        #print words
        m1 = words[1]
        se1 = words[2].split("\n")
        se2 = se1[0]
        #print m1
        #print se2
        meanList.append(m1)
        meanList.append(se2)
    elif words[0] == "22":
        #print words
        m2 = words[1]
        se3 = words[2].split("\n")
        se4 = se3[0]
        #print m2
        #print se4
        meanList.append(m2)
        meanList.append(se4)
    elif words[0] != "":
        if words[0].split("\n")[0] != "g":
            n1 = words[0]
            n2 = words[-2]
            #print words
            #print n1
            #print n2
            meanList.append(n1)
            meanList.append(n2)


# convert list into list of lists for converting into a dataframe

newList = []

i=0
while i < len(meanList):
    newList.append([meanList[i], meanList[i+1], meanList[i+2], meanList[i+3], meanList[i+4], meanList[i+5], meanList[i+6], meanList[i+7], meanList[i+8], meanList[i+9], meanList[i+10] ])
    #print newList
    i = i+11


#create dataframe from list
meanDF = pd.DataFrame.from_records(newList) 


#rename columns
meanDF.columns=["Trait","Chr","Peak_Pos","Peak_Marker1","Peak_Marker2","Mean11","SE11","Mean22","SE22","Num11","Num22"]
meanDF


#write to file
#meanDF.to_csv(path+"/200805_RQTL/jitter/GenoPhenoMeanDF_qtlCon_200828.csv", sep = ",", encoding='utf-8', index=False)


#change numbers to floats to perform calculations
meanDF["Num11"] = meanDF["Num11"].astype(float)
meanDF["Num22"] = meanDF["Num22"].astype(float)
meanDF["Mean11"] = meanDF["Mean11"].astype(float)
meanDF["Mean22"] = meanDF["Mean22"].astype(float)


#perform calculations
meanDF["TotalNum"] = meanDF["Num11"] + meanDF["Num22"]
meanDF["Diff11-22"] = meanDF["Mean11"] - meanDF["Mean22"]
meanDF


# merge files
meanDF = pd.merge(meanDF,parent,on="Trait", how = "outer")
meanDF


# calculate RHE
meanDF["RHE_11-22"] = meanDF["Diff11-22"]/meanDF["Difference"]
meanDF


#write to file
meanDF.to_csv(path+"/200805_RQTL/GenoPhenoMean_RHE_qtlCon_final.csv", sep = ",", encoding='utf-8', index=False)


# file of 11/22 "identity"
genotype = pd.read_csv(path+"/200805_RQTL/jitter/Final_HK_markersList_200828.txt", sep="\t",engine='python') 
genotype.head()


# merge files
mergeGeno = meanDF.merge(genotype, left_on = "Peak_Marker1",right_on = "PeakColq", how = "left")
mergeGeno

# write file
mergeGeno.to_csv(path+"/200805_RQTL/jitter/GenoPhenoMean_RHE_Geno_qtlCon_200828_test.csv", sep = ",", encoding='utf-8')




###########################
# All QTLs (LOCO, R/qtl2) #
###########################


#extract information from text file of the n, mean phenotypic value at peak, etc.
#before running code, double check text files 
#needed to modify to depending on the number outputs and remove the space (if 2 digits instead of 3 digits)
#and whether 22 was placed first

meanList = []
f = open("200817_RQTL2/GenoPhenoMeansNum_10000KLO_200820m.txt")
for line in f:
    words = line.split(" ")
    #print words
    if words[0] == '[1]':
        if words[1][0] == '"':
            key = words[1].split("\n")
            key2 = key[0]
            #print key2
            meanList.append(key2)
        else:
            value = words[1].split("\n")
            #print value[0]
            meanList.append(value[0])
    elif words[0] == "11":
        #print words
        m1 = words[1]
        se1 = words[2].split("\n")
        se2 = se1[0]
        #print m1
        #print se2
        meanList.append(m1)
        meanList.append(se2)
    elif words[0] == "22":
        #print words
        m2 = words[1]
        se3 = words[2].split("\n")
        se4 = se3[0]
        #print m2
        #print se4
        meanList.append(m2)
        meanList.append(se4)
    elif words[0] != "":
        if words[0].split("\n")[0] != "g":
            n1 = words[0]
            n2 = words[-2]
            #print words
            #print n1
            #print n2
            meanList.append(n1)
            meanList.append(n2)
            

#create new list of lists
newList_L = []

i=0
while i < len(meanList):
    newList_L.append([meanList[i], meanList[i+1], meanList[i+2], meanList[i+3], meanList[i+4], meanList[i+5], meanList[i+6], meanList[i+7], meanList[i+8]])
    #print newList
    i = i+9


#create dataframe from list
meanDF_L = pd.DataFrame.from_records(newList_L) 


#rename columns
meanDF_L.columns=["Trait","Chr","Peak_Pos","Mean11","SE11","Mean22","SE22","Num11","Num22"]
meanDF_L


#write to file
#meanDF.to_csv(path+"/200817_RQTL2/GenoPhenoMeanDF_200820.csv", sep = ",", encoding='utf-8', index=False)


#read in file with "marker" information at the peak
peak_marker = pd.read_csv(path+"/200817_RQTL2/peaks_ChromKLO_operm10000_AddEff_200820.txt", sep="\t",engine='python', usecols=[7]) 
peak_marker


# merge peak_marker with rest of dataframe
frames=[meanDF_L, peak_marker]
combDF=pd.concat(frames, axis=1)
combDF


#change numbers to floats to perform calculations
combDF["Num11"] = combDF["Num11"].astype(float)
combDF["Num22"] = combDF["Num22"].astype(float)
combDF["Mean11"] = combDF["Mean11"].astype(float)
combDF["Mean22"] = combDF["Mean22"].astype(float)


#make calculations
combDF["TotalNum"] = combDF["Num11"] + combDF["Num22"]
combDF["Diff11-22"] = combDF["Mean11"] - combDF["Mean22"]
combDF


#merge with parent differences dataframe
combDF = pd.merge(combDF,parent,on="Trait", how = "outer")
combDF


#make RHE calculation
combDF["RHE_11-22"] = combDF["Diff11-22"]/combDF["Difference"]
combDF


#write to file
#combDF.to_csv(path+"/200805_RQTL/GenoPhenoMeanLOD_PVE_PDE_200821.csv", sep = ",", encoding='utf-8', index=False)


#adding some genotypes conversion info (what is 11 and 22) from last time
genotype_L = pd.read_csv(path+"/200817_RQTL2/GenotypePeakInfo_all_200821.txt", sep="\t",engine='python') 
genotype_L.head()


#merge files
mergeGeno_L = combDF.merge(genotype_L, on = "peakCol", how = "left")
mergeGeno_L


#write to file
mergeGeno_L.to_csv(path+"/200817_RQTL2/GenoPhenoMean_RHE_Geno_200821_test.csv", sep = ",", encoding='utf-8')
