Files and script for merging phenotypic datasets into 1 file of the average phenotypic value per individual.

Script:
SumStats_filemerge_clean_v1.R

Input files:
F5random3_NectarALL_200206.txt - nectar - from random 3 flowers per individual (for ind with more than 3 flowers phenotyped)
F5random3_FlowerALL_190913.txt - flower - from random 3 flowers per individual (for ind with more than 3 flowers phenotyped)
NectarySize_NectarOnePlus_merge_200428_mod.txt - nectary - nectary size and associated measurements from all nectaries
F5Seed_180507.txt - seed - all seeds measured before planting
IndSpLine.txt - information on individual, species, line (RiL) numbers

Output file:
output file from script was manually edited (rename column names and remove extra columns)  
becomes the final file to use in calculating correlations, heritability, etc. 

F5_JoinAll3R_IndMean_final_v1.txt - F5 only (removed C & L)
F5_JoinAll3R_LineMeanCL_final_v1.txt - F5 line means (from different output when calculating correlations) and C & L
