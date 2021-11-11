Files and scripts for measuring and calculating nectary size from images.

Dropbox link to nectary photos: https://www.dropbox.com/sh/yuecpvp671tm501/AAAis_xMYj9fz0PymRSoclSza?dl=0
Dropbox link to zipped files: https://www.dropbox.com/s/suo784mu6wtwviy/RILsNectaryPhotos.zip?dl=0
Dryad: https://doi.org/10.5061/dryad.xpnvx0kgz
Raw image files in JPG. Files with "points" designated from FIJI - TIFF (open with FIJI to see the points).

Combine nectary size raw data points from FIJI, and then estimation of nectary size from data points.

Scripts to convert x,y coordinates of points into distances to estimate nectary size:
NectaryDataConversion_for1file_final.py
NectaryDataConversion_mostImages_final.py

Input files:
measurements_191203_200121.txt - use with NectaryDataConversion_for1file_final.py
measurements_19XXXX_20XXXX.txt - use with NectaryDataConversion_mostImages_final.py

Output files:
Calculations_Nectary_allImages_CORRECT_200428.csv
Calculations_Nectary_2-14Images_CORRECT_200428.csv

Script to merge nectar volume and sugar concentration info (and individual info - which nectary is which) to nectary size measurements
Merge_Combine_Nectary_Nectar_final.py

Input files:
F5_RILs_Nectary_Photos.txt - photo number, nectar volume and sugar concentration info
Calculations_Nectary_allImages_CORRECT_200428.csv - oupt put from above script

Output:
Took file from script output, then manually added nectar data to individuals that had 
(1) more than 1 image taken
(2) first 13 photos

File found in Merge_FlowerNectarSeed directory
NectarySize_NectarOnePlus_merge_200428_mod.txt
