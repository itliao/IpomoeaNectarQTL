work_dir <-"~/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles"
setwd (work_dir)

library(qtl)

###############
# genome scan #
###############

#read in data/cross file
Map <- read.cross("csv", work_dir, "RQTL_formatted_200804.csv", estimate.map=FALSE, genotypes = c("11","22")) 

#because data are RILs by selfing
Map <- convert2riself(Map)
summary(Map)
#jitter because many markers with same cM
MapJitter <- jittermap(Map)
summary(MapJitter)

#calculate genotype probabilities
Map <- calc.genoprob(Map, step=1, error.prob=0.001)
MapJitter <- calc.genoprob(MapJitter, step=1, error.prob=0.001)
#genome scan using Haley-Knott, scanone
out.hkJ <- scanone(MapJitter, pheno.col=2:16, method = "hk")

#################################################################
# scantwo permuations - take a while, so done first, then saved #
#################################################################

#set seeds for permutations, split into 2 because it is very time consuming
set.seed(85842518) 
operm2a <- scantwo(mapScan2, pheno.col=2:16, n.perm=500, method="hk") 
save(operm2a, file="perm2eJ.RData")

set.seed(85842519) 
operm2b <- scantwo(mapScan2, pheno.col=2:16, n.perm=500, method="hk")
save(operm2b, file="perm2fJ.RData")

# load output from scantwo permutations 
load("perm2eJ.RData")
load("perm2fJ.RData")
operm2.hkJ <- c(operm2eJ, operm2fJ)

#penalties from 2D permutations
pen <- calc.penalties(operm2.hkJ, alpha = 0.05)

############################################################################
# makeqtl from R/qtl2 LOCO, 1 QTL/chromosome, genome-wide significant QTLs #
# then stepwiseqtl and fitqtl                                              #
############################################################################
#makeqtl from RQTL2 LOCO genome-wide sig 1 peak 
#for each trait
qtlCW <- makeqtl(MapJitter, chr = c(1,2,4,5,8,10,11,13),	pos = c(205,49.931,62.3755,167.7815,206.8125,177.9935,111.0425,65.978), what="prob")
qtlCT <- makeqtl(MapJitter, chr = c(1,7,8,9,13),	pos = c(242.2405,211.835,227.0865,131.927,42.771), what="prob")
qtlCL <- makeqtl(MapJitter, chr = c(1,4,5,8,9,14,15),	pos = c(207,52,169,195.216,145.207,138.9475,27.217), what="prob")
qtlOCL <- makeqtl(MapJitter, chr = c(1,2,4,5,8,9,13,14),	pos = c(207,50,52,169,199.001,145.207,42.771,138.9475), what="prob")
qtlSpL <- makeqtl(MapJitter, chr = c(1,6,7,9,12,14),	pos = c(206,201.194,128.1675,38.0725,144.152,180.5415), what="prob")
qtlLSL <- makeqtl(MapJitter, chr = chr = c(1,4,5,8,9,11,13,14),	pos = c(229.008,46.0925,167.7815,201.2425,145.207,27,58.3455,141.24), what="prob")
qtlPL <- makeqtl(MapJitter, chr = c(1,2,3,4,5,7,8,9,10,13),	pos = c(224.249,112.3145,64.0185,52,167.7815,192.457,215,12.701,58.9175,142.9585), what="prob")
qtlSM <- makeqtl(MapJitter, chr = c(1,12,13),	pos = c(193.618,133.2835,65), what="prob") #actually, no QTLs detected at this level, so using 1 peak, 0.05 chrom
qtlSL <- makeqtl(MapJitter, chr = c(2,12,13),	pos = c(92.3455,116,24.9815), what="prob")
qtlSW <- makeqtl(MapJitter, chr = c(1),	pos = c(193.618), what="prob")
qtlNV <- makeqtl(MapJitter, chr = c(4,5,13),	pos = c(154.002,26,69), what="prob")
qtlNSC <- makeqtl(MapJitter, chr = c(4,9,11,15),	pos = c(197,136.3695,25.8895,60.1785), what="prob")
qtlNS <- makeqtl(MapJitter, chr = c(1,3,4,5,7,8,13,14),	pos = c(204,51,186,52.4275,205.627,30.5875,58.3455,18), what="prob")

#stepwiseqtl conservative
out.sqCWc <- stepwiseqtl(MapJitter, pheno.col=2, max.qtl=11, penalties=pen[1,1:3], qtl=qtlCW, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqCTc <- stepwiseqtl(MapJitter, pheno.col=3, max.qtl=12, penalties=pen[2,1:3], qtl=qtlCT, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqCLc <- stepwiseqtl(MapJitter, pheno.col=4, max.qtl=10, penalties=pen[3,1:3], qtl=qtlCL, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqOCLc <- stepwiseqtl(MapJitter, pheno.col=5, max.qtl=10, penalties=pen[4,1:3], qtl=qtlOCL, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqSpLc <- stepwiseqtl(MapJitter, pheno.col=6, max.qtl=10, penalties=pen[5,1:3], qtl=qtlSpL, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqLSLc <- stepwiseqtl(MapJitter, pheno.col=7, max.qtl=10, penalties=pen[6,1:3], qtl=qtlLSL, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqPLc <- stepwiseqtl(MapJitter, pheno.col=8, max.qtl=12, penalties=pen[7,1:3], qtl=qtlPL, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqSMc <- stepwiseqtl(MapJitter, pheno.col=9, max.qtl=3, penalties=pen[8,1:3], qtl=qtlSM, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqSLc <- stepwiseqtl(MapJitter, pheno.col=10, max.qtl=7, penalties=pen[9,1:3], qtl=qtlSL, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqSWc <- stepwiseqtl(MapJitter, pheno.col=11, max.qtl=2, penalties=pen[10,1:3], qtl=qtlSW, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqNVc <- stepwiseqtl(MapJitter, pheno.col=12, max.qtl=8, penalties=pen[11,1:3], qtl=qtlNV, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqNSCc <- stepwiseqtl(MapJitter, pheno.col=13, max.qtl=10, penalties=pen[12,1:3], qtl=qtlNSC, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)
out.sqNSc <- stepwiseqtl(MapJitter, pheno.col=14, max.qtl=11, penalties=pen[13,1:3], qtl=qtlNS, additive.only = TRUE, method="hk", verbose=2, keeplodprofile = TRUE, keeptrace = TRUE)

#lists for conservative, moderate, liberal stepwiseqtl output
qtlCon <- list(out.sqCWc,out.sqCTc,out.sqCLc,out.sqOCLc,out.sqSpLc,out.sqLSLc,out.sqPLc,out.sqSMc,out.sqSLc,out.sqSWc,out.sqNVc,out.sqNSCc,out.sqNSc)

sink("200805_RQTL/jitter/LOCO1-05g_stepwiseqtl_jitter_seed_200824.txt")
for (i in qtlCon){
  print("StepwiseConAdd (0.05)")
  print(i)
}
sink()

pdf("200805_RQTL/jitter/LOCO1-05g_stepwiseLOD_jitter_seed_200824.pdf")
for(i in qtlCon){
  plotLodProfile(i) 
}
dev.off()

pdf("200805_RQTL/jitter/LOCO1-05g_stepwise_scanone_seed_0.05_200824.pdf")
plotLodProfile(out.sqCWc)
plot(out.hkJ, lodcolumn=1, chr=c(1,2,4,5,8,11,13,14), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqCTc)
plot(out.hkJ, lodcolumn=2, chr=c(1,3,7,8,9,13), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqCLc)
plot(out.hkJ, lodcolumn=3, chr=c(1,4,5,8,9,14,15), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqOCLc)
plot(out.hkJ, lodcolumn=4, chr=c(1,4,5,8,9), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqSpLc)
plot(out.hkJ, lodcolumn=5, chr=c(1,3,5,6,9,12,13,14), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqLSLc)
plot(out.hkJ, lodcolumn=6, chr=c(1,4,5,9,11,13,14), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqPLc)
plot(out.hkJ, lodcolumn=7, chr=c(1,2,3,4,5,7,8,9,10,13), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqSMc)
plot(out.hkJ, lodcolumn=8, chr=c(12), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqSLc)
plot(out.hkJ, lodcolumn=9, chr=c(12,13), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqSWc)
plot(out.hkJ, lodcolumn=10, chr=c(1), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqNVc)
plot(out.hkJ, lodcolumn=11, chr=c(13), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqNSCc)
plot(out.hkJ, lodcolumn=12, chr=c(4,7,9,10,11,15), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqNSc)
plot(out.hkJ, lodcolumn=13, chr=c(1,3,4,5,6,7,8,13,14), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqNVsc)
plot(out.hkJ, lodcolumn=14, chr=c(13), col="red", add=TRUE, lty=2)
plotLodProfile(out.sqNSlc)
plot(out.hkJ, lodcolumn=15, chr=c(1,3,4,5,7,8,13,14), col="red", add=TRUE, lty=2)
dev.off()

#fitqtl summary
sink("200805_RQTL/jitter/LOCO1-05g_fitqtl_seed_0.05_200824.txt")
out.fq <- fitqtl(MapJitter, pheno.col=2, qtl=out.sqCWc, method="hk")
print("CorollaWidth")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=3, qtl=out.sqCTc, method="hk")
print("CorollaThroat")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=4, qtl=out.sqCLc, method="hk")
print("CorollaLength")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=5, qtl=out.sqOCLc, method="hk")
print("OverallCorollaLength")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=6, qtl=out.sqSpLc, method="hk")
print("SepalLength")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=7, qtl=out.sqLSLc, method="hk")
print("LongestStamenLength")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=8, qtl=out.sqPLc, method="hk")
print("PistilLength")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=9, qtl=out.sqSMc, method="hk")
print("SeedMass")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=10, qtl=out.sqSLc, method="hk")
print("SeedLength")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=11, qtl=out.sqSWc, method="hk")
print("SeedWidth")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=12, qtl=out.sqNVc, method="hk")
print("NectarVolume")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=13, qtl=out.sqNSCc, method="hk")
print("SugarConcentration")
print(summary(out.fq))
out.fq <- fitqtl(MapJitter, pheno.col=14, qtl=out.sqNSc, method="hk")
print("NectarySize")
print(summary(out.fq))
sink()

###################################################################
# makeqtl for conservative penalties from output from stepwiseqtl #
###################################################################
#makeqtl for conservative for each trait
qtlCWc <- makeqtl(MapJitter, chr=c(1,2,4,5,8,11,13), pos=c(202.955,49.931,62.376,167.782,227.006,110.718,65.978), what="prob")
qtlCTc <- makeqtl(MapJitter, chr=c(1,1,6,7,8,8,9,13), pos=c(55.517,241.916,60.189,212,192,241.1,131.927,65.978), what="prob")
qtlCLc <- makeqtl(MapJitter, chr=c(1,4,4,5,5,8,9,9,14,15), pos=c(208.535,52,199.142,70,171.078,195.542,27.494,142,139.85,27.217), what="prob")
qtlOCLc <- makeqtl(MapJitter, chr=c(1,4,5,5,8,9), pos=c(204,53,74.937,171.078,195.216,143), what="prob")
qtlSpLc <- makeqtl(MapJitter, chr=c(1,1,3,5,5,5,6,7,9,12,13,14), pos=c(199,254,31.732,63.105,196,199,202.829,128.168,38.073,145,44.336,167.025), what="prob")
qtlLSLc <- makeqtl(MapJitter, chr=c(1,1,4,5,9,11,13,14), pos=c(5,229.662,52,167.782,145.207,27,24.982,141.24), what="prob")
qtlPLc <- makeqtl(MapJitter, chr=c(1,2,3,4,5,7,7,8,8,9,9,10,13,13,15)	, pos=c(206,119.048,59,61,171.078,0,142.411,64.226,214,12.701,135.47,58.918,51,151.345,9), what="prob")
qtlSMc <- makeqtl(MapJitter, chr=c(12), pos=c(146.2), what="prob") 
qtlSLc <- makeqtl(MapJitter, chr=c(12,13), pos=c(115.779,65.978), what="prob")
qtlSWc <- makeqtl(MapJitter, chr=c(1), pos=c(192.64), what="prob")
qtlNVc <- makeqtl(MapJitter, chr=c(13), pos=c(58.346), what="prob")
qtlNSCc <- makeqtl(MapJitter, chr=c(4,7,9,10,11,15), pos=c(197.154,84.8,137,72.608,25.89,60.179), what="prob")
qtlNSc <- makeqtl(MapJitter, chr=c(1,3,4,5,6,7,8,13,13,14,14), pos=c(203,51,182.04145,63.10468,167,205.70889,30.58759,0.24401,142.9588,20,165.87673), what="prob")

#refineqtl for conservative for each trait
RqtlCWc <- refineqtl(MapJitter, pheno.col = 2, qtl=qtlCWc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlCTc <- refineqtl(MapJitter, pheno.col = 3, qtl=qtlCTc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlCLc <- refineqtl(MapJitter, pheno.col = 4, qtl=qtlCLc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlOCLc <- refineqtl(MapJitter, pheno.col =5, qtl=qtlOCLc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlSpLc <- refineqtl(MapJitter, pheno.col = 6, qtl=qtlSpLc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlLSLc <- refineqtl(MapJitter, pheno.col = 7, qtl=qtlLSLc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlPLc <- refineqtl(MapJitter, pheno.col = 8, qtl=qtlPLc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlSMc <- refineqtl(MapJitter, pheno.col = 9, qtl=qtlSMc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlSLc <- refineqtl(MapJitter, pheno.col = 10, qtl=qtlSLc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlSWc <- refineqtl(MapJitter, pheno.col = 11, qtl=qtlSWc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlNVc <- refineqtl(MapJitter, pheno.col = 12, qtl=qtlNVc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlNSCc <- refineqtl(MapJitter, pheno.col = 13, qtl=qtlNSCc, method="hk", verbose=TRUE,keeplodprofile = TRUE)
RqtlNSc <- refineqtl(MapJitter, pheno.col = 14, qtl=qtlNSc, method="hk", verbose=TRUE,keeplodprofile = TRUE)

qtlCon <- list(RqtlCWc,RqtlCTc,RqtlCLc,RqtlOCLc,RqtlSpLc,RqtlLSLc,RqtlPLc,RqtlSMc,RqtlSLc,RqtlSWc,RqtlNVc,RqtlNSCc,RqtlNSc)

sink("200805_RQTL/jitter/Final_refine_jitter_200827.txt")
for (i in qtlCon){
  print("StepwiseConAdd (0.05)")
  print(i)
}
sink()

pdf("200805_RQTL/jitter/Final_refine_jitter_200827.pdf")
for(i in qtlCon){
  plotLodProfile(i) 
}
dev.off()

pdf("200805_RQTL/jitter/Final_refine_scanone_seed_0.05_200826.pdf")
plotLodProfile(RqtlCWc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=1, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlCTc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=2, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlCLc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=3, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlOCLc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=4, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlSpLc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=5, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlLSLc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=6, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlPLc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=7, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlSMc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=8, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlSLc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=9, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlSWc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=10, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlNVc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=11, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlNSCc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=12, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlNSc, showallchr=TRUE)
plot(out.hkJ, lodcolumn=13, col="red", add=TRUE, lty=2)
plotLodProfile(RqtlNVsc, showallchr=TRUE)
dev.off()
sink()

###########################################
# find markers, effect plot, conservative #
###########################################
sink("200805_RQTL/jitter/Final_refine_qtlCon_markersPlots_200826.txt")
pdf("200805_RQTL/jitter/Final_refine_qtlCon_markersPlots_200826.pdf")
a = 1
for (i in qtlCon){
  a <- a+1
  for(j in 1:length(i[[2]])){
    print(j)
    print(a)
    print(i$chr[j])
    print(i$pos[j])
    #Chr <- i$chr[j]
    #Pos <- j$pos[j]
    #print(Chr)
    #print(Pos)
    mar <- find.marker(MapJitter, chr = i$chr[j] , pos = i$pos[j])
    print(mar)
    lodintA <- lodint(i, qtl.index = j, expandtomarkers = TRUE)
    print("LODint_markers")
    print(lodintA)
    bayesintA <- bayesint(i, qtl.index = j, expandtomarkers = TRUE)
    print("BayesInt_markers")
    print(bayesintA)
    effectplot(MapJitter, pheno.col = a, mar )
    plotPXG(MapJitter,pheno.col = a, mar, infer=FALSE)
  }
}
dev.off()
sink()






