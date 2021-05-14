Scripts and files for running QTL analyses in R. 

Scripts:
RQTL_FINAL_stepwise_refineqtl_clean_v1.R - using R/qtl, Haley-Knott regressions
RQTL2_F5_mapJLR_clean_v1.R - using R/qtl2, leave-one-chromosome-out kinship matrices

Input files:
for R/qtl
RQTL_formatted_200804.csv
perm2eJ.RData
perm2fJ.RData

for R/qtl2:
mapJLR_F5_200817.yaml
mapJLR_comp_geno_200615.csv
mapJLR_pheno_200817.csv
mapJLR_gmap_200615.csv
mapJLR_pmap_200615.csv

To convert RQTL2 into RQTL format:
Convert_RQTL2_RQTL_200804_final.py
RQTL_formatted_200804.csv
