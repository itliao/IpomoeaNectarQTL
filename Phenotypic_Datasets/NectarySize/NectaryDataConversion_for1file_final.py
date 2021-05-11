#!/usr/bin/env python
# coding: utf-8


# Goal - to convert the x and y coordinates of points selected on images of nectaries
# to calculate distances between points (e.g. for nectary height, nectary diameter)
# and use these distances to estimate the volume of a nectary modeled after a frustrum

# for just one file (because same image numbers used in other files) 
# see NectaryDataConversion_mostImages_final.py for all the rest

#python conversion from interactive use in jupyter notebook

import pandas as pd
import numpy as np
import random
import re
import os
import math


#path to working directory of file
os.chdir("/Users/irene/Dropbox (Duke Bio_Ea)")
os.chdir("Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/NectaryPhotos-measurement/measurements")


#list files in directory
ls


#read in files
filename = 'measurements_191203_200121.txt'
measure = pd.read_csv(filename, sep="\t", engine='python')


#check file
measure


# only select these 3 columns 
xy = measure[["image","x","y"]]


# check dataframe
xy


#create a dictionary with the key as the image number and the list as the X and Y coordinates of the points
Dict = {}
for row in xy.itertuples():
    key = row[1]
    #print(key)
    if key in Dict:
        #print row[2]
        #print row[3]
        Dict[key].append(row[2])
        Dict[key].append(row[3])
    else:
        Dict[key] = []
        Dict[key].append(row[2])
        Dict[key].append(row[3])
#print Dict


#transform the dictionary into a dataframe
xyDF = pd.DataFrame.from_dict(Dict)
xyDF


#transpose dataframe
xyDFt = pd.DataFrame.transpose(xyDF)
xyDFt


#keep index as a column
xyDFt.reset_index(level=0, inplace=True)
xyDFt


#list to rename columns
ColList=['image', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5', 'x6', 'y6', 'x7', 'y7', 'x8', 'y8']
#rename column names
xyDFt.columns=ColList
xyDFt


#calculate distances between points to get measurements of:
#top outer diameter, top inner diameter, bottom outer diameter
#4 measuremnts of nectary height, 2 measurements of ovary height
xyDFt['Top_out_dia'] = xyDFt['x3'] - xyDFt['x2']
xyDFt['Top_in_dia']= xyDFt['x8'] - xyDFt['x5']
xyDFt['Bottom_out_dia'] = xyDFt['x4'] - xyDFt['x1']
xyDFt['Nec_H1'] = xyDFt['y1'] - xyDFt['y2']
xyDFt['Nec_H2'] = xyDFt['y1'] - xyDFt['y5']
xyDFt['Nec_H3'] = xyDFt['y4'] - xyDFt['y8']
xyDFt['Nec_H4'] = xyDFt['y4'] - xyDFt['y3']
xyDFt['Ov_H1'] = xyDFt['y1'] - xyDFt['y6']
xyDFt['Ov_H2'] = xyDFt['y4'] - xyDFt['y7']


#convert some of the distances to other metrics to use in volume of a frustrum and a cylinder formulas
#average nectary height, top outer radius, top inner radius
#nectary thickness, bottom outer radius, bottom inner radius

#bottom inner radius sub_options:
#A) no thickness. Bottom inner radius == bottom outer radius
#B) 1/2 thickness of the top. Bottom inner radius == bottom outer radiu - 1/2 top nectary thickness
#C) 1/4 thickness of the top. Bottom inner radius == bottom outer radiu - 1/4 top nectary thickness
#D) 1/8 thickness of the top. Bottom inner radius == bottom outer radiu - 1/8 top nectary thickness

#test/calculate all 4, but ended up using option A, bottom inner radius == bottom outer radius

xyDFt['Nec_H_avg'] = (xyDFt['Nec_H1'] + xyDFt['Nec_H2'] + xyDFt['Nec_H3'] + xyDFt['Nec_H4'])/4
xyDFt['Ov_H_avg'] = (xyDFt['Ov_H1'] + xyDFt['Ov_H2'])/2
xyDFt['Top_out_rad'] = xyDFt['Top_out_dia']/2
xyDFt['Top_in_rad'] = xyDFt['Top_in_dia']/2
xyDFt['Nec_thickness'] = xyDFt['Top_out_rad'] - xyDFt['Top_in_rad']
xyDFt['Bottom_out_rad'] = xyDFt['Bottom_out_dia']/2
xyDFt['Bottom_in_rad'] = xyDFt['Bottom_out_rad'] - xyDFt['Nec_thickness']
xyDFt['Bottom_in_rad_A'] = xyDFt['Bottom_out_rad'] 
xyDFt['Bottom_in_rad_B'] = xyDFt['Bottom_out_rad'] - (xyDFt['Nec_thickness']/2)
xyDFt['Bottom_in_rad_C'] = xyDFt['Bottom_out_rad'] - (xyDFt['Nec_thickness']/4)
xyDFt['Bottom_in_rad_D'] = xyDFt['Bottom_out_rad'] - (xyDFt['Nec_thickness']/8)


#calculate each sub-component of calculating the volume of a nectary modeled after the volume of a frustrum or cylinder
xyDFt['FVol_A'] = (math.pi)*(xyDFt['Nec_H_avg'])/3
xyDFt['FVol_B'] = np.power(xyDFt['Top_out_rad'],2) + (xyDFt['Top_out_rad']*xyDFt['Bottom_out_rad']) + np.power(xyDFt['Bottom_out_rad'],2)
xyDFt['FVol_C'] = np.power(xyDFt['Top_in_rad'],2) + (xyDFt['Top_in_rad']*xyDFt['Bottom_in_rad']) + np.power(xyDFt['Bottom_in_rad'],2)
xyDFt['FVol_C1'] = np.power(xyDFt['Top_in_rad'],2) + (xyDFt['Top_in_rad']*xyDFt['Bottom_in_rad_A']) + np.power(xyDFt['Bottom_in_rad_A'],2)
xyDFt['FVol_C2'] = np.power(xyDFt['Top_in_rad'],2) + (xyDFt['Top_in_rad']*xyDFt['Bottom_in_rad_B']) + np.power(xyDFt['Bottom_in_rad_B'],2)
xyDFt['FVol_C3'] = np.power(xyDFt['Top_in_rad'],2) + (xyDFt['Top_in_rad']*xyDFt['Bottom_in_rad_C']) + np.power(xyDFt['Bottom_in_rad_C'],2)
xyDFt['FVol_C4'] = np.power(xyDFt['Top_in_rad'],2) + (xyDFt['Top_in_rad']*xyDFt['Bottom_in_rad_D']) + np.power(xyDFt['Bottom_in_rad_D'],2)
xyDFt['CVol_A'] = (math.pi)*xyDFt['Nec_H_avg']
xyDFt['CVol_B'] = np.power(xyDFt['Top_out_rad'],2) - np.power(xyDFt['Top_in_rad'],2)


#calculate the volume of the nectary
xyDFt['Frustrum_Vol'] = xyDFt['FVol_A']*(xyDFt['FVol_B']-xyDFt['FVol_C'])
xyDFt['Frustrum_Vol_A'] = xyDFt['FVol_A']*(xyDFt['FVol_B']-xyDFt['FVol_C1'])
xyDFt['Frustrum_Vol_B'] = xyDFt['FVol_A']*(xyDFt['FVol_B']-xyDFt['FVol_C2'])
xyDFt['Frustrum_Vol_C'] = xyDFt['FVol_A']*(xyDFt['FVol_B']-xyDFt['FVol_C3'])
xyDFt['Frustrum_Vol_D'] = xyDFt['FVol_A']*(xyDFt['FVol_B']-xyDFt['FVol_C4'])
xyDFt['Cylinder_Vol'] = xyDFt['CVol_A']*xyDFt['CVol_B']


#sort values by image
xyDFt = xyDFt.sort_values(by=['image'])
xyDFt


#write output to file
xyDFt.to_csv('Calculations_Nectary_2-14Images_CORRECT_200428.csv', sep=',', encoding='utf-8', index=False)



