#!/usr/bin/env python
# coding: utf-8

# Goal - to merge nectary size file with nectar data file 
# some potential complication with the first couple of images and individuals that have 
# more than 1 photo, so part of script is to remove those individuals, which are later 
# added back manually

import pandas as pd
import numpy as np
import random
import re
import os
import math


# path to working directory
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/NectaryPhotos-measurement/measurements")


# list files in directory
ls


# file with nectary photo info, nectary info (nectar produced)
nectar = pd.read_csv('F5_RILs_Nectary_Photos.txt', sep = '\t', engine = 'python')
nectar


# only want certain columns
nectar1 = nectar[["Pictures","Ind","Infl_pos","Cor_Len_Var","Nectar_UNC","Sugar_UND","Sugar_D+2.5"]]


# remove the first 12 tester files - interfere with merging
nectar2 = nectar1.iloc[12:2562, :]
nectar2


# create list of individuals with 1+ photos
OnePlus = []
for row in nectar2.itertuples():
    #print row[1] == the photo numbers
    index = row[0]
    if "-" in row[1]:
        OnePlus.append(index)
OnePlus


# dataframe of individuals with 1+ photos
nectarOnePlus = nectar1.iloc[OnePlus, :]
nectarOnePlus


# remove rows with individuals with more than 1 photo
nectar3 = nectar2.drop(OnePlus, axis=0)
nectar3


# rename column names
nectar3 = nectar3.rename(columns ={"Pictures":"image"})
nectar3


# read in file of nectary volume/size calculations
tarySize = pd.read_csv("Calculations_Nectary_allImages_CORRECT_200428.csv")
tarySize


# change "image" column to integers
nectar3['image'] = nectar3['image'].astype(int)


# merge nectar measurements and nectary size measurements
Combo1 = nectar3.merge(tarySize, how="inner", on='image')
Combo1


# write to file
Combo1.to_csv('NectarySize_Nectar_merge_200428.csv', sep=',', encoding='utf-8', index=False)

