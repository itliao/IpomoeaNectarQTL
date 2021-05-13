#!/usr/bin/env python
# coding: utf-8

#python conversion from interactive use in jupyter notebook

# Goal: create files for QTL analyses using rqtl2 
# (1) genetic map
# (2) physical map
# (3) genotypic information, with consensus genotypes between 2 individuals within the same RIL
# (4) for stats - what the percentage of missing genotypes are found for each line

import pandas as pd
import random
import re
import os


# working directory path
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles/FromJoanna")


# list files in directory
ls


# read in geno = genotype file (markers and individuals)
# read in lg = linkage map file (markers, cM, position, chromosome/linkage group)

geno = pd.read_csv('ipomoea_geno_F5_only_200615.csv', sep=",", engine='python')
lg = pd.read_csv('mapJLR_mod_200615.txt', sep="\t",engine='python')


###############################################
# generate genetic and physical map for rqtl2 #
###############################################

#visualize linkage map file
lg.head(10)


# creating a name for the marker - extract the scaffold first
scaf = lg["Scaffold"].str.split("_", n = 2, expand = True)  
lg.insert(2,"Scaf", scaf[1], True)
lg.head(10)


#create the marker name: "scaffold_position"
lg.insert(0, "Marker", lg['Scaf'].map(str) + "_" + lg["Pos"].map(str), True) 
lg.head(10)


# extract chromosome number and add to new column
chrom = lg["Chrom"].str.split("_", n = 2, expand = True)  
lg.insert(1,"Chr", chrom[1], True)
lg.head(10)


#grab the appropriate columns for the genetic map file
gmap1 = lg[["Marker","Chr","cM"]]
gmap1


#write genetic map file for rqtl2
gmap1.to_csv('mapJLR_gmap_200615.csv', sep=",", encoding='utf-8', index=False)


# grab the appropriate columns for the physical map file
pmap = lg[["Marker","Chr","Pos"]]
pmap.head(10)


# write physical map file for rqtl2
pmap.to_csv("mapJLR_pmap_200615.csv", sep=",", encoding='utf-8', index=False)


##################################
# consensus genotypes within RIL #
##################################

#visualize genotype file
geno.head(10)


#create new marker name to match in genetic and physical map
marker = geno["CHR"].str.split("_", n = 3, expand = True)  
geno.insert(0,"Marker", marker[1].map(str)+ "_" + marker[2].map(str), True)
geno.head(10)


# drop "CHR" column with different marker name
genoNew = geno.drop(["CHR"], axis=1)
genoNew.head(10)


# transpose dataframe so that column names == markers and row names == samples
dfT = genoNew.set_index("Marker").transpose()
dfT


# reset index
dfT=dfT.reset_index()
dfT


#get RIL "line" identifiers = the number after "F5_"
sampleName = dfT[["index"]]
sampleName


#split the name to get lines
line = sampleName["index"].str.split("_", n = 3, expand = True) 
line.head(10)


#insert "RIL" line column into the dataframe
dfT.insert(1,"Line", line[2], True)
dfT.head(10)


#drop "index" individual (level_0) column
dfT.drop(dfT.columns[[0]], axis=1, inplace=True)
dfT.head(10)


#sort the dataframe by line
sortedDF = dfT.sort_values(by=['Line'])
sortedDF


#list/set of line names
LineList=[]
LineSet = set()
for i in range(0, len(sortedDF)):
    #LineList.append(testDF.iloc[i,1])
    print(i)
    LineSet.add(sortedDF.iloc[i,0])

for line in LineSet:
    LineList.append(line)

len(LineList) #should expect 322 


#New dictionary with RIL "line" as the key and a list to append to
CombDict = {}
for a in range(len(LineList)):
    key = LineList[a]
    if key in CombDict:
        CombDict[key].append()
    else:
        CombDict[key] = []  

CombDict


#change 21 to be 12 - both reflect heterozygosity
#change from number format to text format in order to do quick-ish comparisons
sortedDF = sortedDF.replace(21,12) #heterozygous
sortedDF2 = sortedDF.replace(12, "12")
sortedDF2 = sortedDF2.replace(11,"11")
sortedDF2 = sortedDF2.replace(22,"22")
sortedDF2.head(10)


#compare genotypes within lines
#likely not the most efficient way

for i in range(len(sortedDF2)-1):
    #print testDF[[i]] <- this prints each column, not row
    #print testDF.iloc[[i]] <- this prints each row
    #print testDF.iloc[[i,2]] <- this prints the first two rows
    #print testDF.iloc[i,1] <- yay! This is what I want - row, column 1 value
    if sortedDF2.iloc[i,0] == sortedDF2.iloc[i+1,0]: #if the lines are equal to each other
        line = sortedDF2.iloc[i,0]
        print(line)
        for j in range(1, len(sortedDF2.columns)):
            if sortedDF2.iloc[i,j] == sortedDF2.iloc[i+1,j]: #if the two genotypes are the same
                geno = sortedDF2.iloc[i,j]
                CombDict[line].append(geno)
            else:
                CombDict[line].append("-")

CombDict


#create empty dataframe with lines as columns
CombDF = pd.DataFrame(columns=tuple(LineList))
print(CombDF)


#check how many genotypes per line - should be 6056, but some are 0
for key, List in CombDict.items():
    #print key
    #print len(List)
    #turns out can't make a dataframe if all index different lengths; need to add singletons later
    if len(List) == 6056:
        print(len(List))
        CombDF[key] = List
    else:
        print(len(List))

CombDF


#make a separate dictionary for RILs that do not have duplicates mainly because
#cannot make a dataframe if all index different lengths; will add these later

SingList=[]
for key, List in CombDict.items():
    #print key
    #print len(List)
    if len(List) != 6056:
        SingList.append(key)

SingList

len(SingList) #should be 28


#combine the marker name with the dataframe with the consensus calls
frame2 = [genoNew[["Marker"]], CombDF]
CombNew = pd.concat(frame2, axis=1)
CombNew


# transpose to switch column and row names
CombNewT = CombNew.transpose()
CombNewT = CombNewT.reset_index()
CombNewT


#drop any rows with NA
CombNewT = CombNewT.dropna(axis=0)
CombNewT


#create dataframe for RILs with only one individual
SingDF = pd.DataFrame()
for i in range(len(sortedDF)):
    for x in SingList:
        if sortedDF2.iloc[i,0] == x:
            print(x)
            AddRow = sortedDF2.iloc[[i]]
            SingDF=SingDF.append(AddRow)

SingDF.head(10)

# recreate columns for combining dataframes downstream 
# not sure why this is necessary, but it works
SingDF.columns = range(SingDF.shape[1])
SingDF

CombNewT.columns = range(CombNewT.shape[1])
CombNewT


# combine dataframe with consensus genotypes of RILs with 2 individuals with genotypes of RILs with 1 individual
frame3 = [CombNewT, SingDF]
TotalDF = pd.concat(frame3)
TotalDF


#rename one of the column names back to marker names, then sort by "Marker" (technically Line)
TotalDF.columns = TotalDF.loc[0]
TotalDF=TotalDF.sort_values(by=["Marker"])
TotalDF = TotalDF.drop([0])
TotalDF

####################################################
# calculate percentage of missing data within line #
####################################################
missDict2 ={}
for i in range(len(LineList)-1):
    miss = TotalDF.iloc[i].value_counts(normalize=True)
    missP = miss.get(key = "-")
    mkey = TotalDF.iloc[i,0]
    missDict2[mkey] = missP
missDict2


missDF2 = pd.DataFrame.from_dict(missDict2, orient='index')
missDF2=missDF2.reset_index()
missDF2 = missDF2.sort_values("index")
missDF2


# write data to files, especially genotype information for rqtl/rqtl2
missDF2.to_csv('missData_JLR_200615.csv', sep=",", encoding='utf-8', index=False)
TotalDF.to_csv('mapJLR_comp_geno_200615.csv', sep=",", encoding='utf-8', index=False)


