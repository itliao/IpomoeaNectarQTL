#!/usr/bin/env python
# coding: utf-8

# python conversion from interactive use in jupyter notebook

# Goal: 
# (1) calculate QTL overlap for two QTL sets (GWS and all QTLs)
# (2) calculating QTL overlap averages within and between modules for both sets

from __future__ import division #need this in order to do division with more fig after decimal point
import pandas as pd
import numpy as np
import math
import random
import re
import os


# path to working directory
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles/")


# get path 
path = os.getcwd()
path


# list files in directory
ls


#read in peaks file
#peakG = QTL peaks from GWS (conservative set) of QTLs merged with output from rqtl and rqtl2
#peakC = QTL peaks from CWS (less conservative set) of QTLs merged with output from rqtl and rqtl2
peakG = pd.read_csv("200817_RQTL2/GenoPheno_Merge_LOCO1g_qtlCon_200902_final.txt", sep="\t",engine='python') 
peakC = pd.read_csv("200817_RQTL2/GenoPhenoMeanLOD_PVE_PDE_Geno_2c05_200901_final.csv", sep=",",engine='python') 


#view files
peakG.head()
peakC.head()


##################
# (1) peakG file #
##################
#create dictionary for each traits
traitD ={}
for index, row in peakG.iterrows():
    trait = row[1]
    if trait not in traitD:
        traitD[trait] = []
#append to lists in dictionary        
for index, row in peakG.iterrows():
    trait = row[1]
    if trait in traitD:
        traitD[trait].append(1)


# examine dictionary of traits (should equal number of QTLs for the trait)
traitD


#create dictionary for lod int for the number of peak overlaps
peakD1 = {}
for indexA, rowA in peakG.iterrows():
    for indexB, rowB in peakG.iterrows():
            traitA = rowA[1]
            traitB = rowB[1]
            key = traitA + "-" + traitB
            if key not in peakD1:
                peakD1[key]=[]
                
peakD1


#iterating through the peakG dataframe
#for LODint (peakD1)
#add 1 if there is an overlap, add 0 if there is no overlap between the two
#For Merge file
for indexA, rowA in peakG.iterrows():
    for indexB, rowB in peakG.iterrows():
            traitA = rowA[1]
            traitB = rowB[1]
            chA = rowA[3]
            chB = rowB[3]
            loA = rowA[6]
            loB = rowB[6]
            hiA = rowA[7]
            hiB = rowB[7]
            key = traitA + "-" + traitB
            print(traitA, traitB, key)
            print(chA, chB)
            print(loA, loB, hiA, hiB)
            if chA == chB:
                if loB >= loA and loB < hiA:
                    peakD1[key].append(1)
                elif  hiB > loA and hiB <= hiA:
                    peakD1[key].append(1)
                elif loB < loA and hiB > hiA:
                    peakD1[key].append(1)
                else:
                    peakD1[key].append(0)


#test sum function for a list within the dictionary
sum(peakD1["CorollaLength-SugarConcentration"])


# dictionary for storing overlap proportion/percentage
overlapD1 = {}
for key in peakD1.keys():
    if key not in overlapD1:
        overlapD1[key] = []
overlapD1


#calculating overlap percentages for LOD int (1)
for key in traitD.keys():
    for key2 in traitD.keys():
        doubleT = key + "-" + key2
        QTL1 = float(sum(traitD[key])) #num of QTL in trait1, as a float
        QTL2 = float(sum(traitD[key2])) #num of QTL in trait2, as a float
        overlap = float(sum(peakD1[doubleT])) #num of overlaping QTLs between trait1 and trait2, as a float
        pOL1 = overlap/QTL1 #percent overlap for trait1
        pOL2 = overlap/QTL2 #percent overlap for trait2
        print(pOL1, pOL2)
        if pOL1 > 1.0:
            pOL1 = 1.0
        if pOL2 > 1.0:
            pOL2 = 1.0
        print(pOL1, pOL2)
        pOverlap = float((pOL1+pOL2)/2) #average overlap between trait1 and trait2
        pOverlap2 = float((overlap)/(QTL1+QTL2-overlap))
        #print(key, key2, doubleT)
        #print(QTL1, QTL2, overlap)
        #print(pOL1, pOL2, pOverlap)
        if doubleT in overlapD1:
            overlapD1[doubleT].append(QTL1)
            overlapD1[doubleT].append(QTL2)
            overlapD1[doubleT].append(overlap)
            overlapD1[doubleT].append(pOL1)
            overlapD1[doubleT].append(pOL2)
            overlapD1[doubleT].append(pOverlap)
            overlapD1[doubleT].append(pOverlap2)
            
overlapD1


# create dataframe from dictionary
overlapDF1 = pd.DataFrame.from_dict(overlapD1)


# transpose dataframe
overlapDF_T1 = overlapDF1.transpose()
overlapDF_T1


# rename columns
overlapDF_T1.columns = ["nQTL1_L", "nQTL2_L","nOverlap_L","OverlapQTL1_L","OverlapQTL2_L","perOverlap_L","Jaccard_Overlap"]


# reset index
overlapDF_T1=overlapDF_T1.reset_index()
overlapDF_T1


# write dataframe to file
overlapDF_T1.to_csv(path + "/200817_RQTL2/QTLOverlap_qtlMerge_200909.csv", sep=",", encoding='utf-8', index=False)


##################
# (2) peakC file #
##################

#create dictionary for traits
traitD ={}
for index, row in peakC.iterrows():
    trait = row[1]
    if trait not in traitD:
        traitD[trait] = []
#append to lists in dictionary        
for index, row in peakC.iterrows():
    trait = row[1]
    if trait in traitD:
        traitD[trait].append(1)

traitD


#create dictionary for lod int for the number of peak overlaps
peakD1 = {}
for indexA, rowA in peakC.iterrows():
    for indexB, rowB in peakC.iterrows():
            traitA = rowA[1]
            traitB = rowB[1]
            key = traitA + "-" + traitB
            if key not in peakD1:
                peakD1[key]=[]

peakD1


#iterating through the peakC dataframe
#for LODint (peakD1)
#add 1 if there is an overlap, add 0 if there is no overlap between the two
#for LOCO file
for indexA, rowA in peakC.iterrows():
    for indexB, rowB in peakC.iterrows():
            traitA = rowA[1]
            traitB = rowB[1]
            chA = rowA[2]
            chB = rowB[2]
            loA = rowA[5]
            loB = rowB[5]
            hiA = rowA[6]
            hiB = rowB[6]
            key = traitA + "-" + traitB
            #print(traitA, traitB, key)
            #print(chA, chB)
            #print(loA, loB, hiA, hiB)
            if chA == chB:
                if loB >= loA and loB < hiA:
                    peakD1[key].append(1)
                elif  hiB > loA and hiB <= hiA:
                    peakD1[key].append(1)
                elif loB < loA and hiB > hiA:
                    peakD1[key].append(1)
                else:
                    peakD1[key].append(0)

peakD1


#test sum function for a pairwise trait list in the dictionary
sum(peakD1["CorollaLength-SugarConcentration"])


# dictionary for storing overlap proportion/percentage
overlapD1C = {}
for key in peakD1.keys():
    if key not in overlapD1C:
        overlapD1C[key] = []


overlapD1C


#calculating overlap percentages for LOD int (1)
for key in traitD.keys():
    for key2 in traitD.keys():
        doubleT = key + "-" + key2
        QTL1 = float(sum(traitD[key])) #num of QTL in trait1, as a float
        QTL2 = float(sum(traitD[key2])) #num of QTL in trait2, as a float
        overlap = float(sum(peakD1[doubleT])) #num of overlaping QTLs between trait1 and trait2, as a float
        pOL1 = overlap/QTL1 #percent overlap for trait1
        pOL2 = overlap/QTL2 #percent overlap for trait2
        print(pOL1, pOL2)
        if pOL1 > 1.0:
            pOL1 = 1.0
        if pOL2 > 1.0:
            pOL2 = 1.0
        print(pOL1, pOL2)
        pOverlap = float((pOL1+pOL2)/2) #average overlap between trait1 and trait2
        pOverlap2 = float((overlap)/(QTL1+QTL2-overlap))
        #print(key, key2, doubleT)
        #print(QTL1, QTL2, overlap)
        #print(pOL1, pOL2, pOverlap)
        if doubleT in overlapD1C:
            overlapD1C[doubleT].append(QTL1)
            overlapD1C[doubleT].append(QTL2)
            overlapD1C[doubleT].append(overlap)
            overlapD1C[doubleT].append(pOL1)
            overlapD1C[doubleT].append(pOL2)
            overlapD1C[doubleT].append(pOverlap)
            overlapD1C[doubleT].append(pOverlap2)


overlapD1C


# create dataframe from dictionary (overlap of CWS QTLs)
overlapDF1C = pd.DataFrame.from_dict(overlapD1C)


# transpose dataframe
overlapDF_T1C = overlapDF1C.transpose()
overlapDF_T1C


# rename columns
overlapDF_T1C.columns = ["nQTL1_L", "nQTL2_L","nOverlap_L","OverlapQTL1_L","OverlapQTL2_L","perOverlap_L","Jaccard_Overlap"]


# reset index
overlapDF_T1C=overlapDF_T1C.reset_index()


# write to file
overlapDF_T1C.to_csv(path + "/200817_RQTL2/QTLOverlap_LOCO2c_200909.csv", sep=",", encoding='utf-8', index=False)


################################################
# calculate average QTL overlap within modules #
################################################

# read in files for pairwise categories 
# noTF indicates no transformed values
pairwiseNoTF = pd.read_csv(path+"/200817_RQTL2/QTLPairwiseCategory_noTF_200902.txt", sep="\t", engine="python")
pairwiseNoTF.head()


#double check overlap files
overlapDF_T1.head() #GWS QTLs
overlapDF_T1C.head() #All QTLs


# for GWS QTLs overlap
mergePWoverlapNoTF_G = overlapDF_T1.merge(pairwiseNoTF, left_on = "index", right_on = "PairwiseTraits", how = "right")
mergePWoverlapNoTF_G


# for CWS QTLs overlap
mergePWoverlapNoTF_C = overlapDF_T1C.merge(pairwiseNoTF, left_on = "index", right_on = "PairwiseTraits", how = "right")
mergePWoverlapNoTF_C


#file in any "NA" with "0.0"
mergePWoverlapNoTF_G=mergePWoverlapNoTF_G.fillna(0.0)
mergePWoverlapNoTF_C=mergePWoverlapNoTF_C.fillna(0.0)


# calculate stats within modules - GWS set
# _L == average percent overlap
# _J == Jaccard Index
meanNoTF_L = mergePWoverlapNoTF_G.groupby("Category")["perOverlap_L"].mean()
stdNoTF_L = mergePWoverlapNoTF_G.groupby("Category")["perOverlap_L"].std()
semNoTF_L = mergePWoverlapNoTF_G.groupby("Category")["perOverlap_L"].sem()

meanNoTF_J = mergePWoverlapNoTF_G.groupby("Category")["Jaccard_Overlap"].mean()
stdNoTF_J = mergePWoverlapNoTF_G.groupby("Category")["Jaccard_Overlap"].std()
semNoTF_J = mergePWoverlapNoTF_G.groupby("Category")["Jaccard_Overlap"].sem()


#convert to frame
meanNoTF_L.df = meanNoTF_L.to_frame()
stdNoTF_L.df = stdNoTF_L.to_frame()
semNoTF_L.df = semNoTF_L.to_frame()

meanNoTF_J.df = meanNoTF_J.to_frame()
stdNoTF_J.df = stdNoTF_J.to_frame()
semNoTF_J.df = semNoTF_J.to_frame()


meanNoTF_J.df


# merge frames together
frames = [meanNoTF_L.df,stdNoTF_L.df,semNoTF_L.df,meanNoTF_J.df,stdNoTF_J.df,semNoTF_J.df]
merge = pd.concat(frames,axis=1)
merge


# rename columns
merge.columns=["meanNoTF_L","stdNoTF_L","semNoTF_L","meanNoTF_J","stdNoTF_J","semNoTF_J"]
merge

# write to file
merge.to_csv(path + "/200817_RQTL2/QTLOverlapSummary_qtlMerge_final.csv", sep=",", encoding='utf-8')


# calculate stats within modules - All/CWS QTLs set
# _L == average percent overlap
# _J == Jaccard Index
meanNoTF_L = mergePWoverlapNoTF_C.groupby("Category")["perOverlap_L"].mean()
stdNoTF_L = mergePWoverlapNoTF_C.groupby("Category")["perOverlap_L"].std()
semNoTF_L = mergePWoverlapNoTF_C.groupby("Category")["perOverlap_L"].sem()

meanNoTF_J = mergePWoverlapNoTF_C.groupby("Category")["Jaccard_Overlap"].mean()
stdNoTF_J = mergePWoverlapNoTF_C.groupby("Category")["Jaccard_Overlap"].std()
semNoTF_J = mergePWoverlapNoTF_C.groupby("Category")["Jaccard_Overlap"].sem()


#convert to frame
meanNoTF_L.df = meanNoTF_L.to_frame()
stdNoTF_L.df = stdNoTF_L.to_frame()
semNoTF_L.df = semNoTF_L.to_frame()

meanNoTF_J.df = meanNoTF_J.to_frame()
stdNoTF_J.df = stdNoTF_J.to_frame()
semNoTF_J.df = semNoTF_J.to_frame()


meanNoTF_J.df


# merge frames together
frames = [meanNoTF_L.df,stdNoTF_L.df,semNoTF_L.df,meanNoTF_J.df,stdNoTF_J.df,semNoTF_J.df]
merge = pd.concat(frames,axis=1)
merge


# rename columns
merge.columns=["meanNoTF_L","stdNoTF_L","semNoTF_L","meanNoTF_J","stdNoTF_J","semNoTF_J"]
merge


# write to file
merge.to_csv(path + "/200817_RQTL2/QTLOverlapSummary_LOCO2g_final.csv", sep=",", encoding='utf-8')

