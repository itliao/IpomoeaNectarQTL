work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/F5_RILs_data/DataAnalysis/CleanedFiles/"
setwd (work_dir)

library(ggplot2)

# plotting two approaches for calculating variance and covariance

VarCovar <- "3R_final3_200803/VarCovar_RILQTL_210929.txt"
lnVar <- read.table(VarCovar, sep="\t", header=TRUE, stringsAsFactor=FALSE)

# plot
p1 <- ggplot(lnVar, aes(x=RIL_LN, y=QTL_LN)) +
  geom_point() + theme_bw() 

pdf("lnVarCovar_210929.pdf")
p1
dev.off()
