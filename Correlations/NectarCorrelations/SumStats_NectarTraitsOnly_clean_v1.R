work_dir <-"/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Lab/F5_RILs/F5_RILs_data/DataAnalysis/CleanedFiles/"
setwd (work_dir)

library(corrplot) #for correlation heatmaps
library(GGally) #for plots
library(car) #for Manova()
library(ggpubr) #for arranging subplots into figure
library(ggpmisc) #for regressions

path <- work_dir

##############
# input file #
##############
file3R <- "NectarTraits_JoinAll3R_IndMean_201110.txt"
file3R <- paste0(path, "3R_final3_200803/", file3R)

F5 <- read.table(file3R, sep="\t", header=TRUE, stringsAsFactor=FALSE)
F5line <- F5$Line

###########################
# Phenotypic correlations #
###########################
#Phenotypic correlations
pcorr <- F5[,-c(1:3)]

#scatterplot of phenotypic correlations
pdf("3R_final3_200803/pCorr_IndMeans3_NectarTraits_Scatterplot_210514.pdf", useDingbats = FALSE)
S2A <- ggpairs(pcorr, 
        upper = list(continuous = wrap("cor", size=2.5)),
        lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
S2A
dev.off()

M <- cor(pcorr)
pdf("3R_final3_200803/pCorr_IndMeans3_210514.pdf", useDingbats = FALSE)
corrplot(M, order = "hclust", method = "color", addrect = 3, hclust.method = "complete", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
corrplot(M, order = "hclust", method = "color", addrect = 3, hclust.method = "ward.D2", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
corrplot(M, order = "hclust", method = "color", addrect = 3, hclust.method = "mcquitty", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
dev.off()

###################################################################
# calculate line means and deviations for correlations (method 1) #
###################################################################
#calculate line means (to calculate genetic correlation)
aggF5 <- aggregate.data.frame(F5, by=list(F5$Line), FUN=mean, na.rm=TRUE)

# merge line means file with individuals file 
F5Ind_Line_Means <- merge(F5, aggF5, by.x="Line", by.y="Group.1")

# find the difference (deviation) between each individual and the line means
# for calculating environmental correlations
# one will be + and one will be -
# first check if the columns to make the calculations are the correct matches
for (i in 4:7){
  print (names(F5Ind_Line_Means)[i])
  print (names(F5Ind_Line_Means)[i+7])
}

# make new columns with deviations
for (i in c(4:7)){
  x <- F5Ind_Line_Means[[i+7]]-F5Ind_Line_Means[[i]]
  #print(x)
  #y <- paste(names(F5Ind_Line_Means[i]),".d",sep="")
  #print(y)
  F5Ind_Line_Means[[i+11]] <- x
}

# create new file for deviations only
F5L_Means_Dev <- F5Ind_Line_Means[,-c(4:10)]

#list of new column names
listNew <-c()
for(i in c(4:7)){
  listNew <- c(listNew,paste(names(F5Ind_Line_Means)[i],".d",sep=""))
}
listNew

#rename deviation columns
names(F5L_Means_Dev)[8:11] <- listNew

#remove rows/lines that there are no deviations (only one line) 
listRm <-c()
for(i in 1:635){
  if (F5L_Means_Dev$NectarVolume.x.d[i] == 0){
    listRm <- c(listRm,i)
  }
}
listRm

#new dataset with just line means 
#313 individuals
New_LMD_rmd <- eF5L_Means_Dev[!duplicated(eF5L_Means_Dev$Line),]

###############################################
# calculate genetic correlations - line means #
###############################################
#genetic correlations
gcorr1 <- New_LMD_rmd[,c(4:7)]

#scatterplot of genetic correlations
pdf("3R_final3_200803/Corr_IndMeans3_Nectar_Scatterplot_210514.pdf", useDingbats = FALSE)
S2B <- ggpairs(gcorr1, 
               upper = list(continuous = wrap("cor", size=2.5)),
               lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
S2B
ggpairs(ecorr1, 
        upper = list(continuous = wrap("cor", size=2.5)),
        lower = list(continuous = wrap("points", alpha = 0.3, size=0.2))) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
dev.off()


############################################################################
# calculating correlations with variance-covariance components from MANOVA #
############################################################################

#variables for main text
fit <- lm(cbind(NectarVolume,SugarConcentation,TotalSugAmt,NectarySize) ~ F5line, data=F5)
manova1 <- Manova(fit)

summan1 <- summary(manova1, multivariate=TRUE)
SSCH <- summan1$multivariate.tests$F5line$SSPH
SSCE <- summan1$multivariate.tests$F5line$SSPE

#write.table(SSCH, file = "3R_final3_200803/Nectar_SSCH_main_3R_201110.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
#write.table(SSCE, file = "3R_final3_200803/Nectar_SSCE_main_3R_201110.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# covariances

#genetic degrees of freedom: 322-1 = 321
#MCPs - line
gMCP1 <- SSCH/321

# environment degrees of freedom == 313
#MCPs - error
eMCP1 <- SSCE/313

gCov1 <- gMCP1 - eMCP1

#convert covariances into correlations
gDF1n <- as.data.frame(cov2cor(gCov1))
eDF1n <- as.data.frame(cov2cor(eMCP1))

#write.table(gDF1n, file = "3R_final3_200803/gCorr_Nectar_SSCH_main_3R_201110.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
#write.table(eDF1n, file = "3R_final3_200803/eCorr_Nectar_SSCE_main_3R_201110.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# genetic correlation 
pdf("3R_final3_200803/gCorr_Nectar_SSCH_3R_main_210514.pdf", useDingbats = FALSE)
M <- as.matrix(gDF1n)
corrplot(M, method = "color", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200)) #plot matrix
corrplot(M, order = "hclust", method="color", addrect = 3, hclust.method = "complete", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "ward.D2", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
corrplot(M, order = "hclust", addrect = 3, method = "color", hclust.method = "mcquitty", col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
corrplot.mixed(M, lower = "color", upper = "number", order = "hclust",hclust.method = "mcquitty", 
               lower.col=colorRampPalette(c("#e6550d", "white", "#660066"))(200),
               upper.col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
S2C <- corrplot.mixed(M, lower = "color", upper = "number", order = "hclust",hclust.method = "mcquitty", 
                      lower.col=colorRampPalette(c("#e6550d", "white", "#660066"))(200),
                      upper.col=colorRampPalette(c("#e6550d", "white", "#660066"))(200))
S2C
dev.off()

############################
# nectar trait regressions #
############################

# all RILs, including RILs with only 1 individual
New_LMD <- F5L_Means_Dev[!duplicated(F5L_Means_Dev$Line),]
write.table(New_LMD, "3R_final3_200803/Nectar_LineMeans_210514.txt", sep = "\t", row.names = FALSE)

# dataset with only relevant variables
ntraits <- New_LMD[,c(4:7)]

lformula <- y ~ x 
qformula <- y ~ x + I(x^2) + 0
# nectary size and nectar volume
S2D <- ggplot(ntraits, aes(NectarySize.y, NectarVolume.y)) +
  geom_point(size=1) +
  stat_smooth(method = "lm", formula = qformula, size = 1, se=FALSE) +
  stat_poly_eq(formula = qformula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse=TRUE) + theme_bw() + labs(x="Nectary Size", y="Nectar Volume") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
S2D

# nectary size and total sugar amount
S2E <- ggplot(ntraits, aes(NectarySize.y, TotalSugAmt.y)) +
  geom_point(size=1) +
  stat_smooth(method = "lm", formula = qformula, size = 1, se=FALSE) +
  stat_poly_eq(formula = qformula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse=TRUE) + theme_bw() + labs(x="Nectary Size", y="Total Sugar Amount") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
S2E

# nectary size and sugar concentration
ts_nv_formula <- (5408-(4442*ntraits$NectarySize.y))/(19.06-(34.29*ntraits$NectarySize.y))
S2Fl <- ggplot(ntraits, aes(NectarySize.y, SugarConcentation.y)) +
  geom_point(size=1) +
  stat_smooth(method = "lm", formula = lformula, size = 1, se=FALSE) +
  stat_poly_eq(formula = lformula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse=TRUE) + theme_bw() + labs(x="Nectary Size", y="Nectar Sugar Concentration") +
  geom_smooth(aes(x=NectarySize.y, y=ts_nv_formula), linetype="dashed", color="red") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

S2Fl

# to obtain coefficients for each model
ntraits$NectarySize2 <- (ntraits$NectarySize.y)^2
nv_q <- lm(NectarVolume.y ~ NectarySize.y + I(NectarySize.y^2) + 0, data=ntraits)
ts_q <- lm(TotalSugAmt.y ~ NectarySize.y + I(NectarySize.y^2) + 0, data=ntraits)
nsc_q <- lm(SugarConcentation.y ~ NectarySize.y + I(NectarySize.y^2) + 0, data=ntraits)

nv_l <- lm(NectarVolume.y ~ NectarySize.y, data=ntraits)
ts_l <- lm(TotalSugAmt.y ~ NectarySize.y, data=ntraits)
nsc_l <- lm(SugarConcentation.y ~ NectarySize.y, data=ntraits)

nv_q
ts_q
nsc_q

nv_l
ts_l
nsc_l

#####################################
# combine regressions into 1 figure #
#####################################

panel1 <- ggarrange(S2D,S2D, ncol=2, labels=c("d)", "d)"))
panel1
panel2 <- ggarrange(S2E,S2Fl, ncol=2, labels=c("e)", "f)"))
panel2

combReg <- ggarrange(panel1, panel2, nrow=2)
combReg

pdf("3R_final3_200803/Nectar_Regressions_210517.pdf", useDingbats = FALSE)
S2D
S2E
S2Fl
combReg
dev.off()




