work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/F5_RILs_data/DataAnalysis/CleanedFiles"
setwd (work_dir)

# merging all phenotypic data files 

#################
# read in files #
#################
# nectar - from random 3 flowers per individual (for ind with more than 3 flowers phenotyped)
nR3all <- read.table("F5random3_NectarALL_200206.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
# flower - from random 3 flowers per individual (for ind with more than 3 flowers phenotyped)
fAll <- read.table("F5random3_FlowerALL_190913.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
# nectary - nectary size and associated measurements from all nectaries
taryAll <- read.table("NectarySize_NectarOnePlus_merge_200428_mod.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
# seed - all seeds measured before planting
seedL <- read.table("F5Seed_180507.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)

#################################
# merge & find individual means #
#################################
# (1) merge files for nectar and nectary
nR3tary <- merge(nR3all, taryAll, by=c("Ind","Cor_Len_Var","Nectar_UNC"))
#write merged tables
#write.table(nR3tary, file = "F5_Nectar(y)_random3all_200803.txt")

# (2) find individual means for nectar and flower
aggn3R <- aggregate.data.frame(nR3tary, by=list(nR3tary$Ind), FUN=mean, na.rm=TRUE)
aggfALL <- aggregate.data.frame(fAll, by=list(fAll$Ind), FUN=mean, na.rm=TRUE)

# (3) re-input information: individual, sp, line numbers
# (3a) read in file with information
indline <- read.table("IndSpLine.txt", sep ="\t", header=TRUE, stringsAsFactor=FALSE)
# (3b) join with individual numbers
n3RL <- merge(aggn3R, indline, by.x="Group.1", by.y="Ind", all.x = TRUE)
fALL_L <- merge(aggfALL, indline, by.x="Group.1", by.y="Ind", all.x = TRUE)

# write output
#write.table(n3RL, file = "F5_Nectar(y)3all_IndMean_200803.txt")
#write.table(fALL_L, file = "F5_FlowerAll_IndMean_200803.txt")

# (4) merge files in following order: 
# flower + seeds
f3RLs <- merge(fALL_L, seedL, by.x="Group.1", by.y="Ind", all.x = TRUE)
# (flower+seeds) + nectar  
All3R <- merge(f3RLs, n3RL, by="Group.1")
# write output
#write.table(All3R, file = "F5_JoinAll3R_IndMean_210510_test.txt")








