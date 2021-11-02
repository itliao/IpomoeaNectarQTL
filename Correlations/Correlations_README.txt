Files and scripts for:
 - calculating line means
 - calculating genetic and environmental correlations using (1) line means and (2) variance-covariance components
 - variance components for calculating heritability
 - comparing the two methods of calculating genetic correlations
 - graphing histograms (Fig. S1)

Script:
SumStats_clean_v2.R - updated 210929 to reflect correction for calculating variances and covariances

Input files:
F5_JoinAll3R_IndMean_final_v1.txt - 635 F5 individuals only
File include the 13 phenotypic values for each individual used in the study.
Flower and nectar traits were averaged from 3 independent measurements. Seed measurements are seed characteristics before planting.

F5_JoinAll3R_LineMeanCL_final_v1.txt - F5 line means and C & L individuals, primarily used for histograms and calculating means and standard errors.
File include the 13 phenotypic values for each RIL and C & L individuals.

CorrCategories_200803.txt
Pairwise traits column name and which module the pairwise traits belong in
The two Cat columns are not used; rather the first column (gCorLinLa) is what is used in the merging of files


Includes a directory for calculating heritability using a python script from output from SumStats_clean_v1.R
