setwd("/Volumes/home/anjohnst/17_member/ATACseq/Quant_summit_2/")


library(stats)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(DESeq2)
library(qvalue)
library(EDASeq)
library(reshape)
source("https://bioconductor.org/biocLite.R")
biocLite("cqn")
library(cqn)
library(scales)
library(reshape2)
options(stringsAsFactors = F, scipen=999)
source("../../Correlation_plots3.R")


count_files <- grep(pattern = "_readPeaks", 
                         x = list.files(), value = TRUE)

rm(list_counts)
list_counts <-list()
df_counts<-NULL
i<-5
head(list_counts[[i]])
for (i in 1:length(count_files)){
  list_counts[[i]] <- read.table(count_files[i], sep="")
  if (i < 2) {
    df_counts <- list_counts[[1]][,]
    #colnames(df_counts)[i+1] <- substr(count_files[i],start = 0, stop = 8)
  }
  else {
    df_counts <- merge(x = df_counts, y = list_counts[[i]][,], by.x = "V2", by.y = "V2", all=T )
  }
}

# there should be 103,476 peaks and 139,588 summits
head(df_counts)
str(df_counts)
dim(df_counts) # 139,588 35

# correctly naming the columns
colnames(df_counts)[-1] <- substr(count_files,start = 0, stop = 9)
class(df_counts[,2])

# replace NAs with 0
sum(is.na(df_counts)) # 33,996
sum(df_counts==0) # NA (b/c no zeros)
df_counts_rmNA <- replace(df_counts, is.na(df_counts), 0)
head(df_counts_rmNA)
sum(df_counts_rmNA==0) #33,996

#generate correlation heatmap for unnormalized data
pdf(file = "ATAC_corr_b4_norm.pdf", width = 8, height = 6)
cormat_b4norm <- generate_heatmap(df_counts_rmNA[,-1],"unNorm peak strength")
dev.off()

sum(rowSums(df_counts_rmNA[,-1])==0) # 0, so there's no peak where every sample is 0.

saveRDS(df_counts_rmNA,"df_counts_rmNA.rds")
# perform quantile normaliztation with GC bias taken into account with CQN

# get the FC content and length (although length is same for all)
GC_df <- 
  read.table("../Total_peak_analysis/Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_GCContent_DF.txt",
             header = T)
head(GC_df)
nrow(GC_df)

# Get the size factor (ie. the # of reads in peaks for each sample)
RiP <- colSums(df_counts_rmNA[,-1])
head(RiP)
nrow(df_counts_rmNA)

cqn.df_counts_rmNA <- cqn(df_counts_rmNA[,-1], lengths = 1,
                  x = GC_df$X8_pct_gc, sizeFactors = RiP,
                  verbose = TRUE, lengthMethod = "fixed")
# Using 'sigma' instead 'sig2' (= sigma^2) is preferred now

saveRDS(cqn.df_counts_rmNA,"cqn.df_counts_rmNA.rds")
#rm(cqn.df_counts_rmNA)

cqn.df_counts_rmNA # 139588
cqnplot(cqn.df_counts_rmNA, n = 1, xlab = "GC content", lty = 1)
names(cqn.df_counts_rmNA)

RPKM.cqn.df_counts <- cqn.df_counts_rmNA$y + cqn.df_counts_rmNA$offset
head(RPKM.cqn.df_counts)

RPM <- sweep(log2(df_counts_rmNA[,-1] + 1), 2, log2(RiP/10^6))
# Warning message:
# In sweep(log2(df_counts_rmNA + 1), 2, log2(RiP/10^6)) :
#  STATS does not recycle exactly across MARGIN
RPKM.std <- sweep(RPM, 1, log2(1 / 10^3))
whGenes <- which(rowMeans(RPKM.std) >= 2 )
length(whGenes)
grp1 <- seq(1,33,2)
grp2 <- seq(2,34,2) 

M.std <- rowMeans(RPKM.std[whGenes, grp1]) - rowMeans(RPKM.std[whGenes, grp2])
A.std <- rowMeans(RPKM.std[whGenes,])
M.cqn <- rowMeans(RPKM.cqn.df_counts[whGenes, grp1]) - rowMeans(RPKM.cqn.df_counts[whGenes, grp2])
A.cqn <- rowMeans(RPKM.cqn.df_counts[whGenes,])

pdf(file = "ATAC_MA_compare_norm.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
       main = "Standard RPKM", ylim = c(-4,4), xlim = c(10,20),
       col = alpha("black", 0.25))
plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
       main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(10,20),
       col = alpha("black", 0.25))
dev.off()  

par(mfrow = c(1,1))
pdf(file = "ATAC_corr_after_CQNnorm.pdf", width = 8, height = 6)
cormat_afterCQN<- generate_heatmap(RPKM.cqn.df_counts,"Norm peak strength")
dev.off()
# only gm83 and 86 don't cluster with their replicate

# filter by the peaks where at least 3 Children have a peak. 

Peak_summit_kids <- read.table("../IDR_peaks/IDR_2017/df_peaks_in_master_kids.txt", sep="\t", header=T)
head(Peak_summit_kids)
dim(Peak_summit_kids) #30363
Peak_summit_inter_summit <- read.table("../Total_peak_analysis/Summit_Summit_peak_inter.txt")
head(Peak_summit_inter_summit)

idx_kids_summits <- which(Peak_summit_inter_summit$V10 %in% Peak_summit_kids$SummitPeak)
length(idx_kids_summits) #57431
summit_kids_names <- Peak_summit_inter_summit[idx_kids_summits,4]
sum(duplicated(summit_kids_names)) #24
head(summit_kids_names)

head(RPKM.cqn.df_counts)
RPKM.cqn.df_counts <- as.data.frame(RPKM.cqn.df_counts)
RPKM.cqn.df_counts$Peak_num <- df_counts_rmNA[,1]
summit_kids_names_num <- do.call(rbind,strsplit(summit_kids_names,"_"))[,2]
head(summit_kids_names_num)
length(summit_kids_names_num)

idx_child_call <- which(RPKM.cqn.df_counts$Peak_num %in% summit_kids_names_num)
length(idx_child_call) # 57,407
head(idx_child_call)

RPKM.cqn.df_counts_child <- RPKM.cqn.df_counts[idx_child_call,]
head(RPKM.cqn.df_counts_child)
dim(RPKM.cqn.df_counts_child) # 57,407    35

pdf(file = "ATAC_corr_after_norm_child.pdf", width = 8, height = 6)
cormat_CQNnorm_child <- generate_heatmap(RPKM.cqn.df_counts_child[,-35],"Norm peak strength child")
dev.off()
# 82 and 86 

## performing the CQN normalization on just the child peaks.
df_counts_child <- df_counts_rmNA[idx_child_call,]
dim(df_counts_child)

GC_df_child <- GC_df[idx_child_call,]
head(GC_df_child)
nrow(GC_df_child)

# Get the size factor (ie. the # of reads in peaks for each sample)
RiP_child <- colSums(df_counts_child[,-1])
head(RiP_child)
length(RiP_child)

cqn.df_child_counts <- cqn(df_counts_child[,-1], lengths = 1,
                          x = GC_df_child$X8_pct_gc, sizeFactors = RiP_child,
                          verbose = TRUE, lengthMethod = "fixed")
saveRDS(cqn.df_child_counts,"cqn.df_child_counts.rds")
rm(cqn.df_child_counts)
rm(A.cqn)
rm(A.std)

cqnplot(cqn.df_child_counts, n = 1, xlab = "GC content", lty = 1)
names(cqn.df_child_counts)

RPKM.cqn.df_child_only_counts<- cqn.df_child_counts$y + cqn.df_child_counts$offset
head(RPKM.cqn.df_child_only_counts)
RPKM.cqn.df_child_only_counts <- as.data.frame(RPKM.cqn.df_child_only_counts)
RPKM.cqn.df_child_only_counts$Peak_num <- df_counts_child[,1]

# perform cluster analysis of data
RPKM.cqn.df_child_only_counts

metadata_merge2 <- metadata_merge
for (i in 1:ncol(metadata_merge)){
  metadata_merge2[,i] <- as.character(metadata_merge[,i])
}
str(metadata_merge2)

Cluster_by_covaritaes <- function(matrix, metadata) {
  d <- dist( t(matrix) )
  library(rafalib)
  mypar()
  hc <- hclust(d)
  
  meta_names <- colnames(metadata)
  for (i in 1:ncol(metadata)){
    plot <- myplclust(hc, labels=colnames(matrix), lab.col=as.fumeric(metadata[,i]),
                      cex=0.5, main=paste("Hierarchical Cluster,\n", meta_names[i],"annotation", sep = " "))
    print(plot)
  }
}
Cluster_by_covaritaes(RPKM.cqn.df_child_only_counts[,-35],metadata_merge2)

str(metadata_merge)


# load in metadata
metadata <- read.table("../Quant_final/metadata.txt",header=T)
head(metadata)

# Before normalization
facts <- as.factor(metadata$sample)
colnames(df_counts_rmNA)
set <- newSeqExpressionSet(as.matrix(df_counts_rmNA[,-1]),
                           phenoData = data.frame(facts, 
                                                  row.names=colnames(df_counts_rmNA)[-1]))

colors <- colorRampPalette(brewer.pal(11, "Spectral"))(17)

sex <- metadata$sex
color_sex <- brewer.pal(3, "Set1")

pdf(file="RLE_ATAC_summit_b4_norm.pdf", width=6, height=4)
plotRLE(set, outline=FALSE, ylim=c(-3, 3), col=colors[facts]) 
title("RLE plot b4 Normalization")
dev.off()

# normalized 
set <- newSeqExpressionSet(as.matrix(RPKM.cqn.df_counts[,-35]),
                           phenoData = data.frame(facts, 
                                                  row.names=colnames(df_counts_rmNA)[-1]))

pdf(file="RLE_ATAC_summit_after_CQN_norm_5x5.pdf", width=6, height=4)
plotRLE(set, outline=FALSE, ylim=c(-.5, .5), col=colors[facts]) 
title("RLE plot after CQN normalization")
dev.off()

pdf(file="PCA_summit_CQN_norm_by_indiv.pdf", width=6, height=6)
plotPCA(set, col=colors[facts], cex=1.2, k=2)
title("PCA plot after CQN normalization by Individual")
dev.off()

pdf(file="PCA_summit_CQN_norm_by_sex.pdf", width=6, height=6)
plotPCA(set, col=color_sex[sex], cex=1.2, k=2)
title("PCA plot after CQN normalization by sex")
dev.off()

# what about normalized subset child peaks 
RPKM.cqn.df_counts_child
set <- newSeqExpressionSet(as.matrix(RPKM.cqn.df_counts_child[,-35]),
                           phenoData = data.frame(facts, 
                                                  row.names=colnames(df_counts_rmNA)[-1]))

pdf(file="RLE_ATAC_summit_after_CQN_norm_5x5_child.pdf", width=6, height=4)
plotRLE(set, outline=FALSE, ylim=c(-.5, .5), col=colors[facts]) 
title("RLE plot after CQN normalization child")
dev.off()

pdf(file="RLE_ATAC_summit_after_CQN_norm_child.pdf", width=6, height=4)
plotRLE(set, outline=FALSE, ylim=c(-3, 3), col=colors[facts]) 
title("RLE plot after CQN normalization child")
dev.off()

pdf(file="PCA_summit_CQN_norm_child_by_indiv.pdf", width=6, height=6)
plotPCA(set, col=colors[facts], cex=1.2, k=2)
title("PCA plot after CQN normalization by Individual")
dev.off()

pdf(file="PCA_summit_CQN_norm_child_by_sex.pdf", width=6, height=6)
plotPCA(set, col=color_sex[sex], cex=1.2, k=2)
title("PCA plot after CQN normalization by sex")
dev.off()

# what about normalized only child peaks 
RPKM.cqn.df_counts_child
RPKM.cqn.df_child_only_counts
set <- newSeqExpressionSet(as.matrix(RPKM.cqn.df_child_only_counts[,-35]),
                           phenoData = data.frame(facts, 
                                                  row.names=colnames(RPKM.cqn.df_child_only_counts)[-35]))

pdf(file="RLE_ATAC_summit_after_CQN_norm_5x5_child_only.pdf", width=6, height=4)
plotRLE(set, outline=FALSE, ylim=c(-.5, .5), col=colors[facts]) 
title("RLE plot after CQN normalization of child only")
dev.off()

pdf(file="RLE_ATAC_summit_after_CQN_norm_child_only.pdf", width=6, height=4)
plotRLE(set, outline=FALSE, ylim=c(-3, 3), col=colors[facts]) 
title("RLE plot after CQN normalization child only")
dev.off()

pdf(file="PCA_summit_CQN_norm_child_only_by_indiv.pdf", width=6, height=6)
plotPCA(set, col=colors[facts], cex=1.2, k=2)
title("PCA plot after CQN normalization by Individual")
dev.off()

pdf(file="PCA_summit_CQN_norm_child_only_by_sex.pdf", width=6, height=6)
plotPCA(set, col=color_sex[sex], cex=1.2, k=2)
title("PCA plot after CQN normalization by sex")
dev.off()

lane1 <- metadata_merge$Lane1
color_lane1 <- brewer.pal(4, "Set1")

pdf(file="PCA_summit_CQN_norm_child_only_by_lane1.pdf", width=6, height=6)
plotPCA(set, col=color_lane1[lane1], cex=1.2, k=2)
title("PCA plot after CQN normalization by lane1")
dev.off()

# rlog transformation? .. Can't because it's already normalized
mat_RPKM_CQN_rlog <- rlog(round(as.matrix(RPKM.cqn.df_counts[,-35])))
# Error in estimateDispersionsFit(object, fitType, quiet = TRUE) : 
# all gene-wise dispersion estimates are within 2 orders of magnitude
# from the minimum value, and so the standard curve fitting techniques will not work.
# One can instead use the gene-wise estimates as final estimates:
#   dds <- estimateDispersionsGeneEst(dds)
# dispersions(dds) <- mcols(dds)$dispGeneEst
# ...then continue with testing using nbinomWaldTest or nbinomLRT

####PCA heatmap plot 
head(metadata)
# adding in the EBV information 
meta2 <- readRDS("../../QTLs_2017/ATAC_QTLs/ATAC_metadata_mean_values.txt")
metadata_merge <- merge(x=metadata, y=meta2[,c(1,8)])
str(meta2)
str(metadata_merge)
metadata_merge$Lane1 <- as.factor(metadata_merge$Lane1)
metadata_merge$Lane2 <- as.factor(metadata_merge$Lane2)
metadata_merge$generation <- as.factor(metadata_merge$generation)
metadata_merge$EBV <- as.factor(metadata_merge$EBV)
metadata_merge$Rep <- as.factor(metadata_merge$Rep)

# get rid of extract []
metadata_merge <- metadata_merge[,-6]
str(metadata_merge)
saveRDS(metadata_merge, "metadata_2.rds")
metadata_merge <- readRDS("../Quant_summit/metadata_2.rds")
####PCA heatmap plot 

x <- princomp(df_counts_rmNA[,-1])
x <- princomp(RPKM.cqn.df_counts[,-35])
x <- princomp(RPKM.cqn.df_counts_child[,-35])
x <- princomp(RPKM.cqn.df_child_only_counts[,-35])

pca <-x 

GROUP.p<-c()
for (i in seq(1:10)){
  GROUP<-metadata_merge[,1]
  GROUP<-lm(pca$loadings[,i]~GROUP)
  GROUP.p<-c(GROUP.p, anova(GROUP)$Pr[1])
}

Replicate.p<-c()
for (i in seq(1:10)){
  Replicate<-metadata_merge[,2]
  Replicate<-lm(pca$loadings[,i]~Replicate)
  Replicate.p<-c(Replicate.p, anova(Replicate)$Pr[1])
}

Sex.p<-c()
for (i in seq(1:10)){
  Sex<-metadata_merge[,3]
  Sex<-lm(pca$loadings[,i]~Sex)
  Sex.p<-c(Sex.p, anova(Sex)$Pr[1])
}

Mito.p<-c()
for (i in seq(1:10)){
  Mito<-metadata_merge[,4]
  Mito<-lm(pca$loadings[,i]~Mito)
  Mito.p<-c(Mito.p, anova(Mito)$Pr[1])
}

harvest.p<-c()
for (i in seq(1:10)){
  harvest<-metadata_merge[,5]
  harvest<-lm(pca$loadings[,i]~harvest)
  harvest.p<-c(harvest.p, anova(harvest)$Pr[1])
}

Lane1.p<-c()
for (i in seq(1:10)){
  Lane1<-metadata_merge[,6]
  Lane1<-lm(pca$loadings[,i]~Lane1)
  Lane1.p<-c(Lane1.p, anova(Lane1)$Pr[1])
}

Lane2.p<-c()
for (i in seq(1:10)){
  Lane2<-metadata_merge[,7]
  Lane2<-lm(pca$loadings[,i]~Lane2)
  Lane2.p<-c(Lane2.p, anova(Lane2)$Pr[1])
}

generation.p<-c()
for (i in seq(1:10)){
  generation<-metadata_merge[,8]
  generation<-lm(pca$loadings[,i]~generation)
  generation.p<-c(generation.p, anova(generation)$Pr[1])
}

EBV.p<-c()
for (i in seq(1:10)){
  EBV<-metadata_merge[,9]
  EBV<-lm(pca$loadings[,i]~EBV)
  EBV.p<-c(EBV.p, anova(EBV)$Pr[1])
}

pvals.raw<-rbind(GROUP.p, Replicate.p,Mito.p,harvest.p,Sex.p,Lane1.p,Lane2.p,generation.p, EBV.p)
pvals.raw<-data.matrix(pvals.raw)
rownames(pvals.raw)<-colnames(metadata_merge)


hmcol3<-colorRampPalette(brewer.pal(9,"Blues"))(5)
hmcol4<-hmcol3
hmcol4[1]<-"white"
logpvals.raw <- -log10(pvals.raw)

# pdf(file="PCA_heatmap_summits_CQN_norm_child.pdf", width=9, height=6)
#pdf(file="PCA_heatmap_summits_b4_norm.pdf", width=9, height=6)
pdf(file="PCA_heatmap_summits_CQN_norm.pdf", width=9, height=6)
heatmap.2(logpvals.raw,Rowv=F,Colv=colnames(pvals.raw),dendrogram='none',trace='none',
          margins=c(8,8),colsep=c(1:11),rowsep=c(1:16),
          sepwidth=c(0.025,0.025), sepcolor="black", col=hmcol4, breaks=c(0,1.30103,2,3,4,max(logpvals.raw)),
          key.xlab = "-log10 pvalue", main="CQN Norm Child ATACseq summit")
dev.off()

## Variance explained by each component:
p.variance.explained <- pca$sdev^2 / sum(pca$sdev^2)
# get top ten
p.variance.explained10 <- p.variance.explained[1:10]
PCs <- rep(0, 10)
for(i in 1:length(PCs)){
  PCs[i] <- as.character(paste("PC ", i, sep = ""))
}
names(p.variance.explained10) <- PCs
# plot percentage of variance explained for each principal component    
#pdf(file="VarBar_summits_CQN_norm_child.pdf", width=6, height=4)
#pdf(file="VarBar_summits_b4_norm.pdf", width=6, height=4)
pdf(file="VarBar_summits_CQN_norm.pdf", width=6, height=4)
barplot(100*p.variance.explained10, las=2, xlab='', ylab='% Variance Explained', 
                main="Variance Explained by component", 
        ylim = c(0,100))
dev.off()


## CQN norm output for QTL analysis

head(RPKM.cqn.df_counts)
dim(RPKM.cqn.df_counts)

RPKM.cqn.df_counts$Peak_num <- paste("Summit_",RPKM.cqn.df_counts$Peak_num,sep="")
head(RPKM.cqn.df_counts$Peak_num)

write.table(RPKM.cqn.df_counts, "Summit_counts_CQN_norm.txt", append = F, quote = F, 
            sep = "\t", row.names = F, col.names = T)

dim(RPKM.cqn.df_counts_child) #57407
head(RPKM.cqn.df_counts_child)
RPKM.cqn.df_counts_child$Peak_num <- paste("Summit_",RPKM.cqn.df_counts_child$Peak_num,sep="")
head(RPKM.cqn.df_counts_child$Peak_num)

write.table(RPKM.cqn.df_counts_child, "Summit_counts_CQN_norm_child.txt", append = F, quote = F, 
            sep = "\t", row.names = F, col.names = T)
saveRDS(RPKM.cqn.df_counts_child, "RPKM.cqn.df_counts_child.rds")

RPKM.cqn.df_counts_child <- readRDS("RPKM.cqn.df_counts_child.rds")
 # Summit_37961 is the IKKB promoter summit
head(RPKM.cqn.df_counts_child)
rownames(RPKM.cqn.df_counts_child) <- RPKM.cqn.df_counts_child$Peak_num
RPKM.cqn.df_counts_child <- RPKM.cqn.df_counts_child[,-35]

library(RFmarkerDetector)
counts_child_rsd <- apply(as.matrix(RPKM.cqn.df_counts_child), 1, rsd)
hist(counts_child_rsd)
boxplot(counts_child_rsd)
summary(counts_child_rsd)
counts_child_rsd[names(counts_child_rsd)=="Summit_37961"]
 # it's in the first quartile
which(names(sort(counts_child_rsd))=="Summit_37961") # 8,657th lowest out of 57407
length(counts_child_rsd) # 57,407

which.min(counts_child_rsd) # is a centromere
 # grep "Summit_82017" Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_merge.bed
hist(as.vector(RPKM.cqn.df_counts_child[RPKM.cqn.df_counts_child$Peak_num=="Summit_37961",-35]))

ATAC_vals <- RPKM.cqn.df_counts_child
head(ATAC_vals)
# need to get the means of ATAC summit values:
ATAC_vals_mean <- data.frame(gm77=(ATAC_vals[,1]+ATAC_vals[,2])/2, 
                             gm78=(ATAC_vals[,3]+ATAC_vals[,4])/2,
                             gm79=(ATAC_vals[,5]+ATAC_vals[,6])/2,
                             gm80=(ATAC_vals[,7]+ATAC_vals[,8])/2,
                             gm81=(ATAC_vals[,9]+ATAC_vals[,10])/2,
                             gm82=(ATAC_vals[,11]+ATAC_vals[,12])/2,
                             gm83=(ATAC_vals[,13]+ATAC_vals[,14])/2,
                             gm84=(ATAC_vals[,15]+ATAC_vals[,16])/2,
                             gm85=(ATAC_vals[,17]+ATAC_vals[,18])/2,
                             gm86=(ATAC_vals[,19]+ATAC_vals[,20])/2,
                             gm87=(ATAC_vals[,21]+ATAC_vals[,22])/2,
                             gm88=(ATAC_vals[,23]+ATAC_vals[,24])/2,
                             gm89=(ATAC_vals[,25]+ATAC_vals[,26])/2,
                             gm90=(ATAC_vals[,27]+ATAC_vals[,28])/2,
                             gm91=(ATAC_vals[,29]+ATAC_vals[,30])/2,
                             gm92=(ATAC_vals[,31]+ATAC_vals[,32])/2,
                             gm93=(ATAC_vals[,33]+ATAC_vals[,34])/2)
rownames(ATAC_vals_mean) <- rownames(ATAC_vals)
head(ATAC_vals_mean)

ATAC_vals_mean

Summit_37961 <- ATAC_vals_mean[rownames(ATAC_vals_mean)=="Summit_37961",c(3:12,17)]
Summit_103041 <- ATAC_vals_mean[rownames(ATAC_vals_mean)=="Summit_103041",c(3:12,17)]

gene_plot_TBC1D4 <- function(x, name) {
  library(ggplot2)
  y <- t(x)
  m <- data.frame(mat= c(1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1))
  p <- data.frame(pat= c(0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0))
  Sex_child_factor <- factor(c("F", "F", "F", "M", "M", "M", "F", "M", "F", "M", "M"))
  temp_df <- cbind(y,m,p, 0,Sex_child_factor)
  colnames(temp_df) <- c("y","m","p","g","sex")
  temp_df <- as.data.frame(temp_df)
  temp_df$Sample <- substr(rownames(temp_df),3,4)
  for (j in 1:nrow(temp_df)) {
    if (temp_df$m[j] == 0 && temp_df$p[j] == 0){temp_df$g[j] <- "0/0"}
    if (temp_df$m[j] == 1 && temp_df$p[j] == 0){temp_df$g[j] <- "1/0"}
    if (temp_df$m[j] == 0 && temp_df$p[j] == 1){temp_df$g[j] <- "0/1"}
    if (temp_df$m[j] == 1 && temp_df$p[j] == 1){temp_df$g[j] <- "1/1"}
  }
  temp_plot <- ggplot(temp_df, aes(factor(g), y, label =Sample))
  temp_plot <- temp_plot + geom_boxplot() + geom_text(aes(color=factor(sex,labels=c("Female","Male")))) +
    geom_jitter(aes(colour=factor(sex,labels=c("Female","Male")))) +
    scale_color_discrete(name = "Gender") +
    ggtitle(name) +
    ylab("Gene Counts") + xlab("Haplotype")
  temp_plot
}

gene_plot_TBC1D4(Summit_37961, "Summit_37961")
gene_plot_TBC1D4(Summit_103041, "Summit_103041")


# how much of EBV amount is an influence on the peak? moderate negative correlation.. 
# but it would make more sense if more EBV led to more NFkB, 
# leading to more signal at the IKKB promoter.
EBV <- substring(meta2$EBV, 2, 5)
plot(EBV[c(3:12,17)], Summit_37961)
cor(x=as.numeric(EBV[c(3:12,17)]), y=as.numeric(Summit_37961))
