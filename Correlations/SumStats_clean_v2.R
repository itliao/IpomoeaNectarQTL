work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/F5_RILs_data/DataAnalysis/CleanedFiles/"
setwd (work_dir)

library(corrplot) #for correlation heatmaps
library(Hmisc) #for rcorr()
library(car) #for qqPlot() and Anova()
library(lme4) #for lmer
library(GGally) #for plots

path <- work_dir

##############
# input file #
##############
# F5 phenotypic information, average mean for each F5 individual 
file3R <- "F5_JoinAll3R_IndMean_final_v1.txt"
F5 <- read.table(file3R, sep="\t", header=TRUE, stringsAsFactor=FALSE)
F5line <- F5$Line

###########################
# Phenotypic correlations #
###########################
# calculated with each F5 individual 
pcorr <- F5[,-c(1:4)]
M <- cor(pcorr)
#pdf("3R_final3_200803/pCorr_IndMeans3_200811.pdf")
corrplot(M, order = "hclust", method = "color", addrect = 3, hclust.method = "complete")
corrplot(M, order = "hclust", method = "color", addrect = 3, hclust.method = "ward.D2")
corrplot(M, order = "hclust", method = "color", addrect = 3, hclust.method = "mcquitty")
#dev.off()
#write.table(M, file = "3R_final3_200803/pCorr_IndMeans3R_200811.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
Res <- rcorr(as.matrix(pcorr), type ="pearson")
#write.table(Res$P, file = "3R_final3_200803/pCorr_IndMeans3R_sig_200803.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
#write.table(Res$r, file = "3R_final3_200803/pCorr_IndMeans3R_200803.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# scatterplot of phenotypic correlations
pdf("3R_final3_200803/pCorr_IndMeans3_Scatterplot_210510_final.pdf", useDingbats=FALSE)
ggpairs(pcorr, 
        upper = list(continuous = wrap("cor", size=2.5)),
        lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
dev.off()

###################################################################
# calculate line means and deviations for correlations (method 1) #
###################################################################

# calculate line means (to calculate genetic correlation)
aggF5 <- aggregate.data.frame(F5, by=list(F5$Line), FUN=mean, na.rm=TRUE)

# merge line means file with individuals file 
F5Ind_Line_Means <- merge(F5, aggF5, by.x="Line", by.y="Group.1")

# find the difference (deviation) between each individual and the line means
# for calculating environmental correlations
# one will be + and one will be -
# first check if the columns to make the calculations are the correct matches
for (i in 5:17){
  print (names(F5Ind_Line_Means)[i])
  print (names(F5Ind_Line_Means)[i+17])
}

# make new columns with deviations
for (i in c(5:17)){
  x <- F5Ind_Line_Means[[i+17]]-F5Ind_Line_Means[[i]]
  #print(x)
  #y <- paste(names(F5Ind_Line_Means[i]),".d",sep="")
  #print(y)
  F5Ind_Line_Means[[i+30]] <- x
}

# create new file for deviations only
F5L_Means_Dev <- F5Ind_Line_Means[,-c(5:21)]

#list of new column names
listNew <-c()
for(i in c(5:17)){
  listNew <- c(listNew,paste(names(F5Ind_Line_Means)[i],".d",sep=""))
}
listNew

# rename deviation columns
names(F5L_Means_Dev)[18:30] <- listNew

# remove rows/lines that there are no deviations (only one line) <- not sure why I would do this..
listRm <-c()
for(i in 1:635){
  if (F5L_Means_Dev$CorollaWidth.x.d[i] == 0){
    listRm <- c(listRm,i)
  }
}
listRm

# new dataset with just duplicates and just deviations for environmental correlations
# 626 individuals
eF5L_Means_Dev <- F5L_Means_Dev[-listRm, ]

# new dataset with of just line means ()
# 313 individuals
New_LMD_rmd <- eF5L_Means_Dev[!duplicated(eF5L_Means_Dev$Line),]

#write.table(New_LMD, file = "3R_final3_200803/R3L_lineMeans_single_200803.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
#write.table(F5L_Means_Dev, file = "3R_final3_200803/R3Lmean_deviations_single_3R_200803.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
#write.table(New_LMD_rmd, file = "3R_final3_200803/R3L_lineMeans_200803.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
#write.table(eF5L_Means_Dev, file = "3R_final3_200803/R3Lmean_deviations_3R_200803.txt", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

#read in new files after adding in the lines that only has one individual in them
#eF5L_Means_Dev <- read.table("3R_final3_200803/R3Lmean_dev_dup_3R_200615m.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
#New_LMD_rmd <- read.table("3R_final3_200803/R3Lmean_dev_3R_200615m.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)

#################################################################
# testing normality with qqPlot and Shapiro-Wilk normality test #
#################################################################
#use line means 
aggF5 <- aggregate.data.frame(F5, by=list(F5$Line), FUN=mean, na.rm=TRUE)

#pdf("3R_final3_200803/qqPlotNormality_200803.pdf")
#sink("3R_final3_200803/NormalityTest_200803.txt")
for (i in c(6:18)){
  x <- colnames(aggF5)[i]
  print(x)
  y <- shapiro.test(aggF5[[i]])
  print(y)
  qqPlot(aggF5[[i]])
}
#sink()
#dev.off()

#####################################################
# calculate genetic correlations - line means       #
# calculate environmental correlations - deviations #
#####################################################
# genetic correlations
gcorr1 <- New_LMD_rmd[,c(5:17)]
M <- cor(gcorr1)
#pdf("3R_final3_200803/gCorr_LineMeans_3R_200811.pdf")
corrplot(M, method = "color") 
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "complete")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "ward.D2")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "mcquitty")
#dev.off()
# write correlations into file
#write.table(M, file = "3R_final3_200803/gCorr_LineMeans_cm_3R_200811.txt")
# another function to calculate correlation and p-value
Res <- rcorr(as.matrix(gcorr1), type ="pearson")
#write.table(Res$P, file = "3R_final3_200803/gCorr_LineMeans_3R_sig_200811.txt")
#write.table(Res$r, file = "3R_final3_200803/gCorr_LineMeans2_3R_200811.txt")

# environmental correlations
ecorr1 <- F5L_Means_Dev[,c(18:30)]
M <- cor(ecorr1)
#pdf("3R_final3_200803/gCorr_LineMeans_3R_200811.pdf")
corrplot(M, method = "color") 
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "complete")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "ward.D2")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "mcquitty")
#dev.off()
# write correlations into file
#write.table(M, file = "3R_final3_200803/eCorr_LineMeans_cm_3R_200811.txt")

# another function to calculate correlation and p-value
Res <- rcorr(as.matrix(ecorr1), type ="pearson")
#write.table(Res$P, file = "3R_final3_200803/eCorr_LineMeans_3R_sig_200811.txt")
#write.table(Res$r, file = "3R_final3_200803/eCorr_LineMeans2_3R_200811.txt")

# cluster dendrograms of correlations #

# genetic correlations
M <- cor(gcorr1)
#pdf("3R_final3_200803/gCorr_LineMeans_3R_cluster_200811.pdf")
d <- as.dist(1-M)
hcorr <- hclust(d, method="complete")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="ward.D2")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="mcquitty")
plot(hcorr)
rect.hclust(hcorr, k=3)
#dev.off()

# environmental correlations
M <- cor(ecorr1)
#pdf("3R_final3_200803/gCorr_LineMeans_3R_cluster_200811.pdf")
d <- as.dist(1-M)
hcorr <- hclust(d, method="complete")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="ward.D2")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="mcquitty")
plot(hcorr)
rect.hclust(hcorr, k=3)
#dev.off()

# scatterplot of genetic and environmental correlations
pdf("3R_final3_200803/Corr_IndMeans3_Scatterplot_210510.pdf", useDingbats=FALSE)
ggpairs(gcorr1, 
        upper = list(continuous = wrap("cor", size=2.5)),
        lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggpairs(ecorr1, 
        upper = list(continuous = wrap("cor", size=2.5)),
        lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
dev.off()

########################################
# calculating broad sense heritability #
########################################
# ANOVA on linear mixed models 
# need to use lmer and random effects to get variance components 

#sink("3R_final3_200803/ANOVA_statsummary_hert_3R_200803.txt")
for (i in c(5:17)){
  trait <- names(F5[i])
  print(trait)
  a <- lm(F5[[i]] ~ F5line, data = F5)
  sumstats <- summary(a)
  aovtest <- Anova(a, type="III")
  print(sumstats)
  print(aovtest)
}
#sink()

# to get variance components - random effect of line using lmer
#sink("3R_final3_200803/ANOVA_hert_variance_3R_200803.txt")
for (i in c(5:17)){
  trait <- names(F5[i])
  print(trait)
  fit.line <- lmer(F5[[i]] ~ (1|F5line), data = F5, REML=TRUE)
  sumstats <- summary(fit.line)
  print(sumstats)
  CI <- confint(fit.line, oldNames = FALSE)
  print(CI)
}
#sink()

############################################################################
# calculating correlations with variance-covariance components from MANOVA #
############################################################################

fit <- lm(cbind(CorollaWidth,CorollaThroat,CorollaLength,OverallCorollaLength,SepalLength,LongestStamenLength,PistilLength,SeedMass,SeedLength,SeedWidth,NectarVolume,SugarConcentration,NectarySize) ~ F5line, data=F5)
manova1 <- Manova(fit)

summan1 <- summary(manova1, multivariate=TRUE)
SSCH <- summan1$multivariate.tests$F5line$SSPH
SSCE <- summan1$multivariate.tests$F5line$SSPE

#write.table(SSCH, file = "3R_final3_200803/SSCH_main_3R_200803.txt")
#write.table(SSCE, file = "3R_final3_200803/SSCE_main_3R_200803.txt")

# covariances

# genetic degrees of freedom: 322-1 = 321
#MCPs - line
gMCP1 <- SSCH/321

# environment degrees of freedom == 313
#MCPs - error
eMCP1 <- SSCE/313

# divide by number of individuals per line, which is 2
gCov1 <- (gMCP1 - eMCP1)/2

write.table(gCov1, file = "3R_final3_200803/gCorr_varCovarUpdate_main_3R_210929.txt")

#convert covariances into correlations
gDF1n <- as.data.frame(cov2cor(gCov1))
eDF1n <- as.data.frame(cov2cor(eMCP1))

# write.table(gDF1n, file = "3R_final3_200803/gCorr_SSCH_main_3R_200803.txt")
# write.table(eDF1n, file = "3R_final3_200803/eCorr_SSCE_main_3R_200803.txt")

# genetic correlation and clusters
M <- gDF1n
#pdf("3R_final3_200803/gCorr_SSCH_3R_main_cluster_200803.pdf")
d <- as.dist(1-M)
hcorr <- hclust(d, method="complete")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="ward.D2")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="mcquitty")
plot(hcorr)
rect.hclust(hcorr, k=3)
M <- as.matrix(gDF1n)
corrplot(M, method = "color") 
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "complete")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "ward.D2")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "mcquitty")
#dev.off()

# environmental correlation and clusters
M <- eDF1n
#pdf("3R_final3_200803/eCorr_SSCE_3R_main_cluster_200803.pdf")
d <- as.dist(1-M)
hcorr <- hclust(d, method="complete")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="ward.D2")
plot(hcorr)
rect.hclust(hcorr, k=3)
hcorr <- hclust(d, method="mcquitty")
plot(hcorr)
rect.hclust(hcorr, k=3)
M <- as.matrix(eDF1n)
corrplot(M, method = "color") 
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "complete")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "ward.D2")
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "mcquitty")
#dev.off()

######################################################
# correlation between the two methods of             #
# calculating genetic and environmental correlations #
######################################################

# convert matrix into 1 column via a list

# genetic correlations - variance-covariance
gDF1L <- c()
for (i in 1:13){
  for (j in 1:13){
    gDF1L <- c(gDF1L, gDF1n[[i]][j])
  }
}
# environmental correlations - variance-covariance
eDF1La <- c()
eDF1Lb <- c()
for (i in 1:13){
  for (j in 1:13){
    eDF1La <- c(eDF1La, paste0(colnames(eDF1n)[i], "-",colnames(eDF1n)[j]))
    eDF1Lb <- c(eDF1Lb, eDF1n[[i]][j])
  }
}
# genetic correlations - line means
gCorLin <- as.data.frame(cor(gcorr1))
gCorLinLa <- c()
gCorLinLb <- c()
for (i in 1:13){
  for (j in 1:13){
    gCorLinLa <- c(gCorLinLa, paste(colnames(gCorLin)[i], "-", colnames(gCorLin)[j]))
    gCorLinLb <- c(gCorLinLb, gCorLin[[i]][j])
  }
}

# environmental correlations - deviations
eCorDev <- as.data.frame(cor(ecorr1))
eCorDevLa <- c()
eCorDevLb <- c()
for (i in 1:13){
  for (j in 1:13){
    eCorDevLa <- c(eCorDevLa, paste(colnames(eCorDev)[i], "-", colnames(eCorDev)[j]))
    eCorDevLb <- c(eCorDevLb, eCorDev[[i]][j])  
  }
}

# genetic covariances 
gCov1 <- c()
for (i in 1:13){
  for (j in 1:13){
    gCov1L <- c(gCov1L, gCov1[[i]][j])
  }
}

# combine two datasets and remove duplicates
CorrCategory <- read.table("CorrCategories_200803.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
#line means vs SSCH
gLM_SS <- cbind(gCorLinLa, gCorLinLb, gDF1L)
gLM_SS <- as.data.frame(gLM_SS)
gLM_SS <- merge(gLM_SS, CorrCategory, by="gCorLinLa")
gLM_SS <- gLM_SS[!duplicated(gLM_SS$gCorLinLb),]
gLM_SS <- gLM_SS[-c(1),]
#write.table(gLM_SS, "3R_final3_200803/gLM_SS_main_200803.txt")

#deviations vs SSCE
eDv_SS <- cbind(gCorLinLa,eCorDevLb, eDF1Lb)
eDv_SS <- as.data.frame(eDv_SS)
eDv_SS <- merge(eDv_SS, CorrCategory, by="gCorLinLa")
eDv_SS <- eDv_SS[!duplicated(eDv_SS$eCorDevLb),]
eDv_SS <- eDv_SS[-c(1),]
#write.table(eDv_SS, "3R_final3_200803/eDV_SS_main_200803.txt")

# plot scatterplots
#pdf("3R_final3_200803/GenEnvCorrPlot_3R_main_200803.pdf", useDingbats=FALSE)
#sink("3R_final3_200803/GenEnvCorr_3R_200803.txt")
gLM_SS$gCorLinLb <- as.numeric(gLM_SS$gCorLinLb) #convert to numeric
gLM_SS$gDF1L <- as.numeric(gLM_SS$gDF1L) #convert to numeric
M1 <- cor(gLM_SS$gCorLinLb,gLM_SS$gDF1L)
print("line means vs SSCH")
print(M1)
p1 <- ggplot(gLM_SS, aes(x=gCorLinLb, y=gDF1L)) +
  geom_point() + theme_bw()
plot(p1)

eDv_SS$eCorDevLb <- as.numeric(eDv_SS$eCorDevLb)
eDv_SS$eDF1Lb <- as.numeric(eDv_SS$eDF1Lb)
M3 <- cor(eDv_SS$eCorDevLb,eDv_SS$eDF1Lb)
print("deviations vs SSCE")
print(M3)
#fig2 <- plot_ly(data = eDv_SSm, x = ~eCorDevLb, y = ~eDF1Lb, color = ~Cat)
p1 <- ggplot(eDv_SS, aes(x=eCorDevLb, y=eDF1Lb)) +
  geom_point() + theme_bw()
plot(p1)

#sink()
#dev.off()

##############
# histograms #
##############
#created a new file with line aggregated (from above)
#also combined with cordatotriloba and lacunosa data
allCLF5 <- "3R_final3_200803/F5_JoinAll3R_LineMeanCL_final_v1.txt"
hisList <- read.table(allCLF5, sep="\t", header=TRUE, stringsAsFactor=FALSE)

#standard error for samples without missing data
ster = function(z) c(se=sd(z) / sqrt((length(z))))
#standard error for samples with missing data (only 1 for some of seed and flowering time traits)
ster2 = function(z) c(se=sd(z, na.rm=TRUE) / sqrt((length(z))-1))

#pdf("3R_final3_200803/F5_3R_histograms_200803.pdf")
#sink("3R_final3_200803/F5_3R_summary_200803.txt")
for (i in c(5:18)){
  x <- names(hisList[i])
  print(x)
  y <- aggregate(hisList[[i]], by=list(hisList$Sp), FUN=mean, na.rm=TRUE)
  a <- aggregate(hisList[[i]], by=list(hisList$Sp), FUN=sd, na.rm=TRUE)
  b <- aggregate(hisList[[i]], by=list(hisList$Sp), FUN=ster)
  c <- aggregate(hisList[[i]], by=list(hisList$Sp), FUN=ster2)
  print(y)
  print(a)
  print(b)
  print(c)
  p <- ggplot(hisList, aes(hisList[[i]], fill=as.factor(Sp), color=as.factor(Sp))) + 
    geom_histogram(position="identity", alpha=0.5) +
    scale_fill_manual(values=c("#CC0066", "#660066", "#333333")) +
    scale_color_manual(values=c("#CC0066", "#660066", "#333333")) +
    geom_vline(data=y, aes(xintercept=x, colour=as.factor(Group.1)), linetype="dashed", size=1) +
    labs(x=names(hisList[i])) + theme_bw() 
  print(p)
}
#sink()
#dev.off()





