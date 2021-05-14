#!/usr/bin/env python
# coding: utf-8

#python conversion from interactive use in jupyter notebook

#Goal: to combine RQTL2 files and necessary info into 1 file to use for R/qtl
#convert "12" to "-" because using riself (RILs), does not take into account heterozygotes 


import pandas as pd
import random
import re
import os


#path to working directory
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles")


# read in files:
# geno = genotypes for each RIL
# pheno = phenotype values for each RIL
# gmap = genetic map (cM) for each marker location
geno = pd.read_csv('mapJLR_comp_geno_200615.csv', sep=",", engine='python')
pheno = pd.read_csv('mapJLR_pheno_200804.txt', sep="\t",engine='python')
gmap = pd.read_csv('mapJLR_gmap_200615.csv', sep=",",engine='python')


# check genotype file
geno.head()


#change 12 to "-" - because heterozygotes are not used for riself QTL analyses
geno2 = geno.replace("12","-")
geno2 = geno2.replace({"Marker":"-"}, "12") #since replace all "12" with "-", need to put the "12" back for the Line
geno2


# check phenotype file
pheno.head()


# merge phenotype and genotype files
geno2Pheno = pd.merge(pheno, geno2, on="Marker", how = "outer")
geno2Pheno


# transpose dataframe
geno2PhenoT = geno2Pheno.transpose().reset_index()
geno2PhenoT


# check genetic map file
gmap.head()


#merge transposed genotype and phenotype file with genetic map file
geno2PhenoTgmap = pd.merge(gmap, geno2PhenoT, left_on="Marker", right_on="index",how="outer")
geno2PhenoTgmap


#attempt to rename columns and drop some extra columns
header = geno2PhenoTgmap['index']
geno2PhenoTgmap.drop(labels=['Marker','index'], axis=1, inplace = True)
geno2PhenoTgmap.insert(0, 'Header', header)
geno2PhenoTgmap


#transpose the dataframe again
geno2PhenoTgmapT = geno2PhenoTgmap.transpose()
geno2PhenoTgmapT


#write to file (still need some manual edits to get to correct format)
geno2PhenoTgmapT.to_csv("RQTL_formatted_200804.csv", sep=",", encoding='utf-8', index=False)






