Files and scripts for: 
(1) calculating QTL overlap using the Jaccard index for both sets of QTLs 
(2) calculating average within and between module QTL overlap.

Script:
QTL_overlap_200909_final.py

Input files:
GenoPheno_Merge_LOCO1g_qtlCon_200902_final.txt
GenoPhenoMeanLOD_PVE_PDE_Geno_2c05_200901_final.csv
QTLPairwiseCategory_noTF_200902.txt

Output files: 
Supporting Information Table S5
Table 2

Some of the QTL 1.5 LOD intervals were modified from the R/qtl or R/qtl2 output to prevent overcounting QTL overlap. 
Modifications are the following:

SepalLength chr 5 peaks 63.1, 196, 199
peak 196, LOD ci_lo: 195 -> 195.89844; Bayes ci_lo: 195.8984 -> 195.89844
peak 199, LOD ci_lo: 196.5509 -> 199.0001; Bayes ci_lo: 196.5509 -> 198.0001

modified ci_lo for
CorollaThroat 8 232.3005 -> 232.3006
PistilLength 4 57.7715 -> 57.7716

