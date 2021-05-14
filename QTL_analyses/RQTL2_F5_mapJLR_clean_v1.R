setwd("~/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles")

library(qtl2)

# Goal - QTL analyses with qtl2, specifically using the LOCO method of identifying QTLs
# using LOCO, one QTL per chromosome, genome-wide significance threshold - part of the GWS QTLs
# using LOCO, multiple QTLs per chromosome, chromosome-wide significance creates the "All QTLs" set

###################
# get files ready #
###################

#create zip
zip_datafiles("mapJLR_F5_200817.yaml", overwrite=TRUE)
#path to zipped files
zipfile <- "mapJLR_F5_200817.zip"
#then read in zip
mapA4 <- read_cross2(zipfile)

##########################################
# follow procedures in R/qtl2 user guide #
##########################################
#step1: calculate genotype probabilities
map <- insert_pseudomarkers(mapA4$gmap, step=1)
pr <- calc_genoprob(mapA4, map, error_prob=0.001)

#leave one chromosome out - LOCO method
#scan chrom using kinship matrix data from all other chrom
kin_loco <- calc_kinship(pr, "loco")

#perform genome scan with loco
out_klo <- scan1(pr, mapA4$pheno, kin_loco)

#output - matrix of LOD scores, positions x phenotypes
#plot() to plot LOD curves
#lodcolumn - column/trait to plot
par(mar=c(5.1, 4.1, 1.1, 1.1))

#plot genome scan results per trait with different y-axis limits
pdf("QTLs_mapJLR_LOCO_201216.pdf")
for (i in 1:13){
  print(i)
  plot(out_klo, map, lodcolumn = i, col="slateblue", ylim=c(0, 5))
  plot(out_klo, map, lodcolumn = i, col="violetred", ylim=c(0, 10))
  plot(out_klo, map, lodcolumn = i, col="darkgreen", ylim=c(0, 20))
}
dev.off()

#perform permutation test - establish stats significant
#n_perm = number of permutation replicates.
operm4 <- scan1perm(pr, mapA4$pheno, kin_loco, n_perm=1000)
sum4 <- summary(operm4, alpha=c(0.2, 0.1, 0.05, 0.01))

save(operm4, file="opermLOCO_RQTL2.RData")
write.table(sum4, file = "200817_RQTL2/sum4perm_mapJLRklo1000_200817.txt")

#find_peaks() - id peaks
#LOD support, use drop 1.5 with peakdrop 1.8
#find the 0.05 genome-wide threshold for each trait
peaksCW <- find_peaks(out_klo, map, threshold=3.258233314, peakdrop=1.8, drop=1.5)
peaksCT <- find_peaks(out_klo, map, threshold=3.271583122, peakdrop=1.8, drop=1.5)
peaksCL <- find_peaks(out_klo, map, threshold=3.229508436, peakdrop=1.8, drop=1.5)
peaksCLV <- find_peaks(out_klo, map, threshold=3.109889231, peakdrop=1.8, drop=1.5)
peaksSL <- find_peaks(out_klo, map, threshold=3.073749999, peakdrop=1.8, drop=1.5)
peaksMFL <- find_peaks(out_klo, map, threshold=3.284543513, peakdrop=1.8, drop=1.5)
peaksPL <- find_peaks(out_klo, map, threshold=3.225302135, peakdrop=1.8, drop=1.5)
peaksSM <- find_peaks(out_klo, map, threshold=3.252094823, peakdrop=1.8, drop=1.5)
peaksSdL <- find_peaks(out_klo, map, threshold=3.108473759, peakdrop=1.8, drop=1.5)
peaksSW <- find_peaks(out_klo, map, threshold=3.180948084, peakdrop=1.8, drop=1.5)
peaksNV <- find_peaks(out_klo, map, threshold=3.223716159, peakdrop=1.8, drop=1.5)
peaksNSC <- find_peaks(out_klo, map, threshold=3.174653448, peakdrop=1.8, drop=1.5)
peaksNS <- find_peaks(out_klo, map, threshold=3.237713343, peakdrop=1.8, drop=1.5)

sink("200817_RQTL2/peaks_LOD_LOCO_operm1000_200817.txt")
print("CorollaWidth")
print(peaksCW)
print("CorollaThroat")
print(peaksCT) 
print("CorollaLength")
print(peaksCL) 
print("OverallCorollaLength")
print(peaksCLV)
print("SepalLength")
print(peaksSL) 
print("LongestStamenLength")
print(peaksMFL) 
print("PistilLength")
print(peaksPL) 
print("SeedMass")
print(peaksSM) 
print("SeedLength")
print(peaksSdL)
print("SeedWidth")
print(peaksSW) 
print("NectarVolume")
print(peaksNV) 
print("SugarConcentration")
print(peaksNSC) 
print("NectarySize")
print(peaksNS)
sink()

###########################################
# chromosome-level significance with LOCO #
###########################################
# create leave-one-chromosome-out kinship matrix 
kinship1 <- calc_kinship(pr, "loco")[[1]]
kinship2 <- calc_kinship(pr, "loco")[[2]]
kinship3 <- calc_kinship(pr, "loco")[[3]]
kinship4 <- calc_kinship(pr, "loco")[[4]]
kinship5 <- calc_kinship(pr, "loco")[[5]]
kinship6 <- calc_kinship(pr, "loco")[[6]]
kinship7 <- calc_kinship(pr, "loco")[[7]]
kinship8 <- calc_kinship(pr, "loco")[[8]]
kinship9 <- calc_kinship(pr, "loco")[[9]]
kinship10 <- calc_kinship(pr, "loco")[[10]]
kinship11 <- calc_kinship(pr, "loco")[[11]]
kinship12 <- calc_kinship(pr, "loco")[[12]]
kinship13 <- calc_kinship(pr, "loco")[[13]]
kinship14 <- calc_kinship(pr, "loco")[[14]]
kinship15 <- calc_kinship(pr, "loco")[[15]]

# find chromosome-level significance for loco with 10,000 permutations
operm1_klo <- scan1perm(pr[,"1"], mapA4$pheno, kinship1, n_perm=10000)
operm2_klo <- scan1perm(pr[,"2"], mapA4$pheno, kinship2, n_perm=10000)
operm3_klo <- scan1perm(pr[,"3"], mapA4$pheno, kinship3, n_perm=10000)
operm4_klo <- scan1perm(pr[,"4"], mapA4$pheno, kinship4, n_perm=10000)
operm5_klo <- scan1perm(pr[,"5"], mapA4$pheno, kinship5, n_perm=10000)
operm6_klo <- scan1perm(pr[,"6"], mapA4$pheno, kinship6, n_perm=10000)
operm7_klo <- scan1perm(pr[,"7"], mapA4$pheno, kinship7, n_perm=10000)
operm8_klo <- scan1perm(pr[,"8"], mapA4$pheno, kinship8, n_perm=10000)
operm9_klo <- scan1perm(pr[,"9"], mapA4$pheno, kinship9, n_perm=10000)
operm10_klo <- scan1perm(pr[,"10"], mapA4$pheno, kinship10, n_perm=10000)
operm11_klo <- scan1perm(pr[,"11"], mapA4$pheno, kinship11, n_perm=10000)
operm12_klo <- scan1perm(pr[,"12"], mapA4$pheno, kinship12, n_perm=10000)
operm13_klo <- scan1perm(pr[,"13"], mapA4$pheno, kinship13, n_perm=10000)

save(operm1_klo, file="opermLOCO_chr1_RQTL2.RData")
save(operm2_klo, file="opermLOCO_chr2_RQTL2.RData")
save(operm3_klo, file="opermLOCO_chr3_RQTL2.RData")
save(operm4_klo, file="opermLOCO_chr4_RQTL2.RData")
save(operm5_klo, file="opermLOCO_chr5_RQTL2.RData")
save(operm6_klo, file="opermLOCO_chr6_RQTL2.RData")
save(operm7_klo, file="opermLOCO_chr7_RQTL2.RData")
save(operm8_klo, file="opermLOCO_chr8_RQTL2.RData")
save(operm9_klo, file="opermLOCO_chr9_RQTL2.RData")
save(operm10_klo, file="opermLOCO_chr10_RQTL2.RData")
save(operm11_klo, file="opermLOCO_chr11_RQTL2.RData")
save(operm12_klo, file="opermLOCO_chr12_RQTL2.RData")
save(operm13_klo, file="opermLOCO_chr13_RQTL2.RData")

sink("200817_RQTL2/sumChromPerm10000_KLO_200819.txt")
summary(operm1_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm2_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm3_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm4_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm5_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm6_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm7_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm8_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm9_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm10_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm11_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm12_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
summary(operm13_klo, alpha=c(0.2, 0.1, 0.05, 0.01))
sink()

#peaks with max and min (0.05) at the chrom level
#multiple peaks on a chromosome
#peakdrop indicates the amount that the LOD curve must drop below the lowest of two adjacent peaks
peaksMaxChrom <- find_peaks(out_klo, map, threshold=2.14, drop=1.5)
peaksMaxChrom2 <- find_peaks(out_klo, map, threshold=2.14, peakdrop=1.8, drop=1.5)
peaksMinChrom <- find_peaks(out_klo, map, threshold=1.73, drop=1.5)
peaksMinChrom2 <- find_peaks(out_klo, map, threshold=1.73, peakdrop=1.8, drop=1.5)

sink("200817_RQTL2/peaks_ChromKLO_operm10000_200819.txt")
print("Max of max, 1 peak")
print(peaksMaxChrom) 
print("Max of max, 2 peak")
print(peaksMaxChrom2)
print("Min of min, 1 peak")
print(peaksMinChrom)
print("Min of min, 2 peak")
print(peaksMinChrom2) 
sink()

############################################################
# to find markers at confidence interval locations (in cM) #
############################################################
#read in files with peaks information

peaksKLO <- read.table("200817_RQTL2/peaks_forMarkers_KLO_all_200820.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
peaksHK <- read.table("200817_RQTL2/peaks_forMarkers_HK_all_200820.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)

#use "map" to find markers - these include the pseudomarkers that were inserted
pmar <- find_marker(map, 4, 74)

#iterate through the dataframe to get the markers
peakList <- list()
loList <- list()
hiList <- list()

for (i in 1:nrow(peaksKLO)){
  chr <- peaksKLO[i,3]
  peak <- peaksKLO[i,4]
  ci_lo <- peaksKLO[i,6]
  ci_hi <- peaksKLO[i,7]
  pmarP <- find_marker(map, chr, peak)
  pmarL <- find_marker(map, chr, ci_lo)
  pmarH <- find_marker(map, chr, ci_hi)
  peakList[[i]] <- pmarP
  loList[[i]] <- pmarL
  hiList[[i]] <- pmarH
}

#make the list into a dataframe column
peakCol <- do.call(rbind,peakList)
loCol <- do.call(rbind,loList)
hiCol <- do.call(rbind,hiList)
#add/append the "datafames" into the main dataframe
peaksKLO <- cbind(peaksKLO,peakCol)
peaksKLO<- cbind(peaksKLO,loCol)
peaksKLO <- cbind(peaksKLO,hiCol)
#write to new table
write.table(peaksKLO, file = "peaks_ChromKLO_operm10000_markers_200820.txt",  append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

# get closest physical location of pseudomarkers 
# (interpolated positions for the pseudomarkers)
tointerp <- list("1" = c(pos1.204=204, pos1.205=205, pos1.206=206, pos1.207=207, pos1.69=69, pos1.70=70),
                 "2" = c(pos2.50=50),
                 "3" = c(pos3.141=141, pos3.32=32, pos3.51=51, pos3.52=52),
                 "4" = c(pos4.186=186, pos4.197=197, pos4.224=224, pos4.52=52, pos4.61=61),
                 "5" = c(pos5.169=169, pos5.202=202, pos5.26=26, pos5.27=27, pos5.98=98),
                 "6" = c(pos6.159=159, pos6.57=57),
                 "7" = c(pos7.43=43),
                 "8" = c(pos8.215=215),
                 "9" = c(pos9.28=28, pos9.35=35),
                 "11" = c(pos11.27=27),
                 "12" = c(pos12.116=116, pos12.119=119, pos12.84=84),
                 "13" = c(pos13.142=142, pos13.6=6, pos13.65=65, pos13.69=69, pos13.26=26),
                 "14" = c(pos14.163=163, pos14.18=18, pos14.182=182, pos14.184=184, pos14.96=96, pos14.97=97),
                 "15" = c(pos15.10=10, pos15.7=7)
                 )

sink("200817_RQTL2/PseudoMarkers_KLO_200820.txt")
interp_map(tointerp, mapA4$gmap, mapA4$pmap)
sink()

###################################################################################
# find means, n, genotype info of QTL peaks from all QTLs, chromosome-level, LOCO #
# use information to calculate RHE                                                #
###################################################################################

#re-read in the file (couldn't figure out how to drop the levels in the original file)
peaksKLOnew <- read.table("200817_RQTL2/peaks_ChromKLO_operm10000_markers_200820.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)

pdf("200817_RQTL2/GenoPheno_10000KLO_200820.pdf")
sink("200817_RQTL2/GenoPhenoValues_10000KLO_200820.txt")
for (i in 1:nrow(peaksKLOnew)){
  chrom <- peaksKLOnew[i,3]
  peak <- peaksKLOnew[i,4]
  trait <- peaksKLOnew[i,2]
  pL <- peaksKLOnew[i,8]
  pLa <- find_marker(mapA4$gmap, chrom, peak)
  g <- maxmarg(pr, map, chr=chrom, pos=peak, return_char = TRUE)
  g2 <- maxmarg(pr, mapA4$gmap, chr=chrom, pos=peak, return_char = TRUE)
  p <- plot_pxg(g, mapA4$pheno[,trait], ylab = trait, SEmult=2, force_labels = TRUE, omit_points = FALSE)
  num <- table(g)
  print(trait)
  print(chrom)
  print(peak)
  print(pL)
  print(g)
  print(pLa)
  print(g2)
}
sink()
dev.off()

sink("200817_RQTL2/GenoPhenoMeansNum_10000KLO_200820.txt")
for (i in 1:nrow(peaksKLOnew)){
  chrom <- peaksKLOnew[i,3]
  peak <- peaksKLOnew[i,4]
  trait <- peaksKLOnew[i,2]
  pL <- peaksKLOnew[i,8]
  pLa <- find_marker(mapA4$gmap, chrom, peak)
  g <- maxmarg(pr, map, chr=chrom, pos=peak, return_char = TRUE)
  g2 <- maxmarg(pr, mapA4$gmap, chr=chrom, pos=peak, return_char = TRUE)
  p <- plot_pxg(g, mapA4$pheno[,trait], ylab = trait, SEmult=2, force_labels = TRUE, omit_points = FALSE)
  num <- table(g)
  print(trait)
  print(chrom)
  print(peak)
  print(pL)
  print(p)
  print(num)
}
sink()

#################################################################################
# find means, n, genotype info of markers from R/qtl final stepwise set of QTLs #
# use the conservative (0.05, HK) set of QTLs                                   #
# use information to calculate RHE                                              #
#################################################################################
peaksHKFinal <- read.csv("200805_RQTL/jitter/Final_qtlCon_markers_pos_200827.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)

#single putative QTL position, plot_pxg()
#vector of genotypes from maxmarg()
#return_char = TRUE, vector of character strings with genotype labels
#creates plots 
pdf("200817_RQTL2/Final_GenoPheno_qtlCon_200827.pdf")
sink("200817_RQTL2/Final_GenoPhenoValues_qtlCon_200827.txt") #only for marker and genotypes, don't print means and numbers
for (i in 1:nrow(peaksHKFinal)){
  chrom <- peaksHKFinal[i,3]
  peak <- peaksHKFinal[i,4]
  trait <- peaksHKFinal[i,2]
  pL <- peaksHKFinal[i,5]
  pLa <- find_marker(mapA4$gmap, chrom, peak)
  g <- maxmarg(pr, map, chr=chrom, pos=peak, return_char = TRUE)
  g2 <- maxmarg(pr, mapA4$gmap, chr=chrom, pos=peak, return_char = TRUE)
  p <- plot_pxg(g, mapA4$pheno[,trait], ylab = trait, SEmult=2, force_labels = TRUE, omit_points = FALSE)
  num <- table(g)
  print(trait)
  print(chrom)
  print(peak)
  print(pL) #peak marker
  print(g) #all the genotypes
  print(pLa) #peak marker
  print(g2) # all the genotypes
}
sink()
dev.off()

sink("200817_RQTL2/Final_GenoPhenoMeansNum_qtlCon_200827.txt") #only for means and numbers, don't print pL and g
for (i in 1:nrow(peaksHKFinal)){
  chrom <- peaksHKFinal[i,3]
  peak <- peaksHKFinal[i,4]
  trait <- peaksHKFinal[i,2]
  pL <- peaksHKFinal[i,5]
  pLa <- find_marker(mapA4$gmap, chrom, peak)
  g <- maxmarg(pr, map, chr=chrom, pos=peak, return_char = TRUE)
  g2 <- maxmarg(pr, mapA4$gmap, chr=chrom, pos=peak, return_char = TRUE)
  p <- plot_pxg(g, mapA4$pheno[,trait], ylab = trait, SEmult=2, force_labels = TRUE, omit_points = FALSE)
  num <- table(g)
  print(trait)
  print(chrom)
  print(pL) #peak marker
  print(peak)
  print(p) #print the means
  print(num) #print the numbers of 11 and 22
}
sink()

###################################################
# create summary "barplots" for both sets of QTLs #
###################################################

#read in the peak files from qtlMerge (GWS QTLs) and LOCO2c (all QTLs)
qtlMerge <- read.table("200817_RQTL2/peaks_qtlMerge_201004.csv", sep=",",header=TRUE, stringsAsFactor=FALSE)
LOCO2c <- read.table("200817_RQTL2/peaks_LOCO2c_201004.csv", sep=",", header=TRUE, stringsAsFactor=FALSE)

pdf("PeaksSummary_201004.pdf")
plot_peaks(qtlMerge, map)
plot_peaks(LOCO2c, map)
dev.off()
