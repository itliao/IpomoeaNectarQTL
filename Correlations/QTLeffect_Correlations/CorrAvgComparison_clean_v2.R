work_dir <- ("/Users/irene/Dropbox (Duke Bio_Ea)/Rausher Lab/Current_lab_members/Irene_Liao_Data/Dissertation/1_QTL_chapter/Dissertation/Mark_Analyses/Correlations/")
setwd(work_dir)

library(GGally) #for plots
library(ggpubr) #to arrange subfigures into figure

# Goal: 
# (1) Find correlations between QTL overlap, genetic correlations, bias, predicted genetic correlations
# (2) Plot correlations

# path to working directory
path <- work_dir

# read in files
# CWS_all - values from all QTL pairwise traits
# CWS_FN - values from only flower and nectar QTLs (all QTLs)
# GWS_all - values from GWS QTL pairwise traits
# GWS_FN - values from only flower and nectar QTLs (GWS QTLs)

CWS_all <- read.table("CWS_all_201005.csv", sep=",", header=TRUE, stringsAsFactor=FALSE)
CWS_FN <- read.table("CWS_FN_201005.csv", sep=",", header=TRUE, stringsAsFactor=FALSE)
GWS_all <- read.table("GWS_all_201005.csv", sep=",", header=TRUE, stringsAsFactor=FALSE)
GWS_FN <- read.table("GWS_FN_201005.csv", sep=",", header=TRUE, stringsAsFactor=FALSE)

# examine the correlations for all values
pdf("Correlations_allTypes_201004.pdf", useDingbats = FALSE)
ggpairs(CWS_all[2:8])
ggpairs(CWS_FN[2:7])
ggpairs(GWS_all[2:8])
ggpairs(GWS_FN[2:8])
dev.off()

###########
# CWS_all #
###########
pdf("Correlations_CWS_all_210514.pdf", useDingbats = FALSE)
#sink("Correlations_CWS_all_201005.txt")
M1 <- cor(CWS_all$gen_corr,CWS_all$pred_corr)
print("(1) gen_corr vs pred_corr")
print(M1)
p1 <- ggplot(CWS_all, aes(x=gen_corr, y=pred_corr)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
A2 <- ggplot(CWS_all, aes(x=gen_corr, y=pred_corr, color=Cat2)) +
  geom_point(size=1) + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(A2)

M1 <- cor(CWS_all$gen_corr,CWS_all$Correct_QTL_overlap)
print("(2) gen_corr vs Correct_QTL_overlap")
print(M1)
B4 <- ggplot(CWS_all, aes(x=gen_corr, y=Correct_QTL_overlap, color=Cat1)) +
  geom_point(size=1) + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap (all QTLs)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(B4)
p1 <- ggplot(CWS_all, aes(x=gen_corr, y=Correct_QTL_overlap, color=Cat2)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)

M1 <- cor(CWS_all$avg_total_RHE,CWS_all$bias)
print("(3) avg total RHE vs bias")
print(M1)
p1 <- ggplot(CWS_all, aes(x=avg_total_RHE, y=bias)) +
  geom_point() + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
C2 <- ggplot(CWS_all, aes(x=avg_total_RHE, y=bias, color=Cat2)) +
  geom_point(size=1) + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(C2)
#sink()
dev.off()

##########
# CWS_FN #
##########
pdf("Correlations_CWS_FN_210514.pdf", useDingbats = FALSE)
#sink("Correlations_CWS_FN_201005.txt")
M1 <- cor(CWS_FN$gen_corr,CWS_FN$pred_corr)
print("(1) gen_corr vs pred_corr")
print(M1)
B2 <- ggplot(CWS_FN, aes(x=gen_corr, y=pred_corr)) +
  geom_point(size=1) + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(B2)
p1 <- ggplot(CWS_FN, aes(x=gen_corr, y=pred_corr, color=Cat1)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)


M1 <- cor(CWS_FN$gen_corr,CWS_FN$Correct_QTL_overlap)
print("(2) gen_corr vs Correct_QTL_overlap")
print(M1)
p1 <- ggplot(CWS_FN, aes(x=gen_corr, y=Correct_QTL_overlap)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
p1 <- ggplot(CWS_FN, aes(x=gen_corr, y=Correct_QTL_overlap, color=Cat1)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)

M1 <- cor(CWS_FN$avg_total_RHE,CWS_FN$bias)
print("(3) avg total RHE vs bias")
print(M1)
D2 <- ggplot(CWS_FN, aes(x=avg_total_RHE, y=bias)) +
  geom_point(size=1) + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(D2)
p1 <- ggplot(CWS_FN, aes(x=avg_total_RHE, y=bias, color=Cat1)) +
  geom_point() + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)

#sink()
dev.off()

###########
# GWS_all #
###########
pdf("Correlations_GWS_all_210514.pdf", useDingbats = FALSE)
#sink("Correlations_GWS_all_201005.txt")
M1 <- cor(GWS_all$gen_corr,GWS_all$correct_pred_corr)
print("(1) gen_corr vs correct_pred_corr")
print(M1)
p1 <- ggplot(GWS_all, aes(x=gen_corr, y=correct_pred_corr)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
A1 <- ggplot(GWS_all, aes(x=gen_corr, y=correct_pred_corr, color=Cat2)) +
  geom_point(size=1) + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(A1)
#fig1

M1 <- cor(GWS_all$gen_corr,GWS_all$Correct_QTL_overlap)
print("(2) gen_corr vs Correct_QTL_overlap")
print(M1)
A4 <- ggplot(GWS_all, aes(x=gen_corr, y=Correct_QTL_overlap, color=Cat1)) +
  geom_point(size=1) + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap (GWS QTLs)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(A4)
p1 <- ggplot(GWS_all, aes(x=gen_corr, y=Correct_QTL_overlap, color=Cat2)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)

M1 <- cor(GWS_all$Correct_ave_total_RHE,GWS_all$correct_bias)
print("(4) correct_ave_total_RHE vs correct_bias")
print(M1)
p1 <- ggplot(GWS_all, aes(x=Correct_ave_total_RHE, y=correct_bias)) +
  geom_point() + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
C1 <- ggplot(GWS_all, aes(x=Correct_ave_total_RHE, y=correct_bias, color=Cat2)) +
  geom_point(size=1) + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(C1)

#sink()
dev.off()

##########
# GWS_FN #
##########
pdf("Correlations_GWS_FN_210514.pdf", useDingbats = FALSE)
#sink("Correlations_GWS_FN_201005.txt")
M1 <- cor(GWS_FN$gen_corr,GWS_FN$correct.pred_corr)
print("(1) gen_corr vs correct_pred_corr")
print(M1)
B1 <- ggplot(GWS_FN, aes(x=gen_corr, y=correct.pred_corr)) +
  geom_point(size=1) + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(B1)
p1 <- ggplot(GWS_FN, aes(x=gen_corr, y=correct.pred_corr, color=Cat1)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="predicted correlations (rQ)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)


M1 <- cor(GWS_FN$gen_corr,GWS_FN$Correct.QTL_overlap)
print("(2) gen_corr vs Correct_QTL_overlap")
print(M1)
p1 <- ggplot(GWS_FN, aes(x=gen_corr, y=Correct.QTL_overlap)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
p1 <- ggplot(GWS_FN, aes(x=gen_corr, y=Correct.QTL_overlap, color=Cat1)) +
  geom_point() + theme_bw() + labs(x="genetic correlations (rG)", y="QTL overlap") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)

M1 <- cor(GWS_FN$Correct.ave_total_RHE,GWS_FN$correct.bias)
print("(3) correct_ave_total_RHE vs correct_bias")
print(M1)
D1 <- ggplot(GWS_FN, aes(x=Correct.ave_total_RHE, y=correct.bias)) +
  geom_point(size=1) + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(D1)
p1 <- ggplot(GWS_FN, aes(x=Correct.ave_total_RHE, y=correct.bias, color=Cat1)) +
  geom_point() + theme_bw() + labs(x="average total RHE", y="bias (rQ - rG)") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot(p1)
#sink()
dev.off()

####################################
# arrange to make bioxriv Figure 4 #
####################################
pdf("Correlations_Figure4_arrange.pdf", useDingbats = FALSE)
fig4 <- ggarrange(A4, B4, ncol = 2, nrow=2, labels=c("a)","b)"))
fig4
dev.off()

pdf("Correlations_genCorrQTLOverlap_moduleColors.pdf", useDingbats = FALSE)
A4
B4
dev.off()

############################
# arrange to make Figure 5 #
############################
pdf("Correlations_Figure5_arrange.pdf", useDingbats = FALSE)
panel1 <- ggarrange(A1, A2, ncol = 2, nrow=1, labels=c("a)","b)"), legend = "none")
panel1
panel2 <- ggarrange(B1, B2, ncol = 2, nrow=1, labels=c("c)","d)"))
panel2
panel3 <- ggarrange(C1, C2, ncol = 2, nrow=1, labels=c("e)","f)"), legend = "none")
panel3
panel4 <- ggarrange(D1, D2, ncol = 2, nrow=1, labels=c("g)","h)"))
panel4

comb <- ggarrange(panel1, panel2, panel3, panel4,
                  nrow=4)
comb
dev.off()