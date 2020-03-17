# Script for data analyses of Keeffe & Blackburn:
#       Load volume files and place landmarks
#       Plot PCA
#       Run ANOVA
#       Plot phylomorphospace
#       Run phylogenetic regressions
#       Plot supplemental figure
#       Plot ASE figure
#
# Fixed landmarks used: 17
#

#Loading needed packages

rm(list=ls(all=T))
library(geomorph)
library(ape)
library(phytools)
library(geiger)
library(calibrate)
library(RevGadgets)
library(dplyr)
setwd("D:/Desktop/Humerus_test")

#Importing volume files and placing landmarks; example file

Adenomus_kelartii <- read.ply("Adenomus_kelartii_1580_right_humerus.ply")
digit.fixed(Adenomus_kelartii, fixed = 17, index=F, ptsize=1, center=T)

#Compiling generated landmark files

filelist <- list.files(pattern = ".nts")
mydata <- readmulti.nts(filelist)

#Running generalized procrustes and PCA on compiled landmark files, extracting PCA scores

Y.gpa <- gpagen(mydata)
PCA <- plotTangentSpace(Y.gpa$coords, axis1 = 1, axis2 = 2, label=TRUE)
Eigen <- PCA$sdev^2				
VarPCs <- round(100*Eigen/sum(Eigen),2)
scores <- as.data.frame(PCA$pc.scores)
x <- as.list(scores$PC1)
y <- as.list(scores$PC2)

#Loading grouping variables

classifier <- read.csv("classifier.csv", header=T, row.names=1)

#Plotting PC scores 1 and 2

plot(x, y,
     xlab = "PC1 (37.64%)",
     ylab = "PC2 (21.67%)",
     main = "Forward Burrower vs. Non-burrower Humeri",
     cex = 2,
     pch = c(15,18,17,16)[classifier$Sex],
     col = c("blue","green","forestgreen","black")[classifier$Type])
abline(v=0, h=0, lty=5)
labs <- filelist
textxy(x,y,labs)

#Re-running PCA on trimmed dataset in preparation for ANOVA 
#       unspecified-burrowing species removed
#       groupings reduced to forward and non-forward burrowing

binarybins_no_unspecified <- read.csv("binarybins_no_unspecified.csv", header=T, row.names=1)
cut_filelist_vector <- as.vector(binarybins_no_unspecified$File)
mydata_cut <- readmulti.nts(cut_filelist_vector)
Y.gpa <- gpagen(mydata_cut)
PCA <- plotTangentSpace(Y.gpa$coords, axis1 = 1, axis2 = 2, label=TRUE)
Eigen <- PCA$sdev^2				
VarPCs <- round(100*Eigen/sum(Eigen),2)
scores <- as.data.frame(PCA$pc.scores)
x <- as.list(scores$PC1)
y <- as.list(scores$PC2)

#Plotting trimmed PCA to check for errors

plot(x, y,
     xlab = "PC1 (37.64%)",
     ylab = "PC2 (21.67%)",
     main = "Forward Burrower vs. Non-burrower Humeri",
     cex = 2,
     pch = c(15,18,17,16)[binarybins_no_unspecified$Sex],
     col = c("forestgreen","black")[binarybins_no_unspecified$Type])
abline(v=0, h=0, lty=5)
labs <- cut_filelist_vector
textxy(x,y,labs)

#Preparing data for ANOVA tests

data(mydata_cut)
Y.gpa <- gpagen(mydata_cut)
gdf <- geomorph.data.frame(Y.gpa, burr_type = binarybins_no_unspecified$Type)

#ANOVA using a simple allometry model

fit.size <- procD.lm(coords ~ log(Csize), data = gdf, iter=10000)
fit.anova <- anova(fit.size)

#ANOVA between burrowing/non-burrowing bins

Burr_Anova <- procD.lm(coords ~ burr_type, data=gdf, iter=10000) 
results.anova <- anova(Burr_Anova)

#Reading Jetz and Pyron 2018 tree for phylomorphospace

frog.tree <- read.nexus("amph_shl_new_Consensus_7238.phy")

#Pruning tips

species_list <- read.csv("phylomorph_pruning.csv", row.names=1)
burr.mode <- setNames(species_list[,1],rownames(species_list))
tips <- rownames(species_list)
namecheck <- name.check(frog.tree, burr.mode, data.names=tips)
treenotdata <- as.vector(namecheck$tree_not_data)
datanottree <- as.vector(namecheck$data_not_tree)
droppedtip <- drop.tip(frog.tree,namecheck$tree_not_data)

#Plotting trimmed tree to check for errors

plotTree(droppedtip,type="fan",fsize=0.7,ftype="i",lwd=1) 
cols <- setNames(c("blue","green","gray","forestgreen"),levels(burr.mode))
tiplabels(pie=to.matrix(burr.mode[droppedtip$tip.label],
                        levels(burr.mode)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

#Preparing tree for phylomorphospace

tt <- multi2di(droppedtip)
tt$edge.length[tt$edge.length==0]<-1e-8

#Preparing data for phylomorphospace

phylomorph_bins <- read.csv("phylomorph_bins.csv", header=T)
phylomorph_bins_vector <- as.vector(phylomorph_bins)
cols <- c("blue","green","gray","forestgreen")[phylomorph_bins_vector$Burrowing.Type]

#Re-running PCA on trimmed dataset in preparation for phylomorphospace
#       copied landmark files to subfolder
#       duplicate species removed
#       landmark files renamed to match tree tips

setwd("D:/Desktop/Humerus_test/phylomorph_landmarks")
filelist <- list.files(pattern = ".nts")
mydata <- readmulti.nts(filelist)
data(mydata)
Y.gpa <- gpagen(mydata)

#Plotting phylomorphospace
dev.off()
plotGMPhyloMorphoSpace(tt,Y.gpa$coords, 
                       node.labels = FALSE,
                       ancStates = FALSE, 
                       plot.param = list(t.bg=cols, txt.cex=0.75))

#Calculating the degree of phylogenetic signal from Procrustes shape variables

physignal(Y.gpa$coords, phy = tt, iter = 10000)

#Plotting Supplimental figure 1

setwd("D:/Desktop/Humerus_test")
frog.tree <- read.nexus("amph_shl_new_Consensus_7238.phy")
genera_bins <- read.csv("frog_genera_bins.csv", row.names=1, header=F)
burr.mode <- setNames(genera_bins[,1],rownames(genera_bins))
tips <- rownames(genera_bins)
namecheck <- name.check(frog.tree, burr.mode, data.names=tips)
treenotdata <- as.vector(namecheck$tree_not_data)
datanottree <- as.vector(namecheck$data_not_tree)
droppedtip_genera <- drop.tip(frog.tree,namecheck$tree_not_data)
genera_bins <- read.csv("frog_genera_bins.csv", header=F)
states_list <- as.factor(genera_bins$V2)
names_vector <- as.vector(genera_bins$V1)
genera_characters <- setNames(states_list, names_vector)
b <- names(genera_characters)[genera_characters=="b"]
c <- names(genera_characters)[genera_characters=="c"]
d <- names(genera_characters)[genera_characters=="d"]
tt <- paintBranches(droppedtip_genera,edge=sapply(b,match,droppedtip_genera$tip.label),
                  state="b",anc.state="a")
tt <- paintBranches(tt,edge=sapply(c,match,droppedtip_genera$tip.label),
                  state="c")
tt <- paintBranches(tt,edge=sapply(d,match,droppedtip_genera$tip.label),
                  state="d")
cols <- setNames(c("black","blue","red", "green"),c("a","b","c","d"))
plot(tt,colors=cols,lwd=1,split.vertical=TRUE,ftype="i",type="fan")

#Plotting ancestral state tree generated by .REV file

freeK_RJ_tree_file = "ancestral_states_ase_freeK_RJ.tree"
freeK_RJ <- plot_ancestral_states(freeK_RJ_tree_file, summary_statistic="MAP",
                                  tip_label_size=0.5,
                                  xlim_visible=NULL,
                                  node_label_size=0,
                                  show_posterior_legend=TRUE,
                                  node_size_range=c(1,3),
                                  alpha=0.75)
plot(freeK_RJ)