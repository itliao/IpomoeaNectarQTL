setwd("~/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Scripts_Instructions/QTL/RQTLfiles")

library(qtl2)

# Goal - to find the closest marker and number of CC and LL genotypes
# at the midpoint location between 2 overlapping QTLs

###################
# get files ready #
###################

#path to zipfiles
zipfile <- "mapJLR_F5_200817.zip"

#read in zip files
mapA4 <- read_cross2(zipfile)

#calculate genotype probabilities
map <- insert_pseudomarkers(mapA4$gmap, step=1)
pr <- calc_genoprob(mapA4, map, error_prob=0.001)

####################################################################################################
# find info: closest marker to midpoint, number of genotypes at marker, identity of 11/22 genotype #
####################################################################################################

#file of midpoint values between 2 overlapping QTLs
#for all QTLs, use file: Midpoints_MDR_200922.csv
midpoint <- read.csv("200817_RQTL2/Midpoints_MDR_200922.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
#midpoint <- read.csv("200817_RQTL2/Midpoints_qtlMerge_200924.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)

midpoint[,5:8]<-NA
for (i in 1:nrow(midpoint)){
  chrom <- midpoint[i,3]
  peak <- midpoint[i,4]
  pLa <- find_marker(mapA4$gmap, chrom, peak)
  markerP <- find_markerpos(mapA4, pLa, na.rm = TRUE)
  g <- maxmarg(pr, map, chr=chrom, pos=peak, return_char = TRUE)
  g2 <- maxmarg(pr, mapA4$gmap, chr=chrom, pos=peak, return_char = TRUE)
  num <- table(g)
  midpoint[i,5]<-pLa
  midpoint[i,6]<-markerP[2][[1]]
  midpoint[i,7]<-num[1][[1]]
  midpoint[i,8]<-num[2][[1]]
}

#rename column names
names(midpoint)[5] <- "Marker"
names(midpoint)[6] <- "marker_genetic_pos"
names(midpoint)[7] <- "Num11"
names(midpoint)[8] <- "Num22"

midpoint$Scaffold <- sapply(strsplit(midpoint$Marker, "_"), head, 1)

# read in file with 11/22 identity by scaffold
geno1122 <- read.csv("200817_RQTL2/Genotypes_Scaffold.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
geno1122 <- geno1122[,2:4]

#merge files
mergedFile <- merge(midpoint, geno1122, by = "Scaffold")

#write to file
#write.table(mergedFile,"200817_RQTL2/Midpoints_qtlMerge.txt",sep="\t", row.names = FALSE)
write.table(mergedFile,"200817_RQTL2/Midpoints_allQTLs.txt",sep="\t", row.names = FALSE)
