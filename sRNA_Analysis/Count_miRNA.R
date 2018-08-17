setwd("/Volumes/home-1/anjohnst/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Novel_miRNA/Count_miRNA_novel/")

library(EDASeq)
library(gplots)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(reshape)
library(ggplot2)
library(edgeR)

options(stringsAsFactors = F)

## Load and merge count files to create master table
miRNA_Files <- grep("miRNA.txt.summary", 
                    list.files(), value =TRUE, invert = TRUE)
miRNA_Files <- miRNA_Files[-c(1:3,21:25)]
miRNA_list <- list()
i<-2
head(miRNA_list)
head(miRNA_list[[i]])
nrow(miRNA_list[[3]])
for (i in 1:length(miRNA_Files)){
  miRNA_list[[i]] <- read.table(miRNA_Files[i], header = T)
  miRNA_list[[i]] <- miRNA_list[[i]][,c(1,2,7)]
  if (i < 2 ) {
    Mat_miRNA_norm <-  miRNA_list[[i]]
    #Mat_miRNA_novel <- miRNA_list[[i]][Mat_miRNA_Novelidx,c(1,3)]
    #Mat_miRNA_known <- miRNA_list[[i]][-Mat_miRNA_Novelidx,c(1,3)]
  }
  else {
    Mat_miRNA_norm <- merge(Mat_miRNA_norm, miRNA_list[[i]][,c(1,3)], 
                            by = "Geneid", all= TRUE) 
    #Mat_miRNA_novel <- merge(Mat_miRNA_novel,miRNA_list[[i]][Mat_miRNA_Novelidx,c(1,2)], 
    #                        by = "Geneid") 
  }
}


rownames(Mat_miRNA_norm)<- Mat_miRNA_norm[,1]
Mat_miRNA_norm <- Mat_miRNA_norm[,-1]
head(Mat_miRNA_norm)
tail(Mat_miRNA_norm)
nrow(Mat_miRNA_norm) #3,034
saveRDS(Mat_miRNA_norm, "Mat_miRNA_norm.rds")

summary(Mat_miRNA_norm[,5])

sum(rowSums(Mat_miRNA_norm[grep("novel", rownames(Mat_miRNA_norm), invert=T),-1])>1)

sum(rowMeans(Mat_miRNA_norm[grep("novel", rownames(Mat_miRNA_norm), invert=T),-1])>5)
# 618

idx_expres <-which(rowMeans(Mat_miRNA_norm[grep("novel", rownames(Mat_miRNA_norm), invert=T),-1])>5)
sum(rowSums(Mat_miRNA_norm[grep("novel", rownames(Mat_miRNA_norm), invert=T),c(4:13,18)]>5)>3) # 632

# Filter out sex chromosome
no_x <- Mat_miRNA_norm[Mat_miRNA_norm$Chr!= "chrX",]
noSex <- no_x[no_x$Chr !="chrY",]
nrow(noSex)
# 2832 vs. 2756
str(noSex)

# Filter out miRNA that are not expressed (<5 counts) in at least 3 of the children
idx_samp <- rowSums(noSex[,c(4:13,18)]>5)>3
length(idx_samp)
sum(idx_samp) #722

miRNA_noXY_samp <- noSex[idx_samp,]
nrow(miRNA_noXY_samp) # 722 vs. 1143
head(miRNA_noXY_samp)
saveRDS(miRNA_noXY_samp, "miRNA_noXY_exp_thres.rds")

unique(miRNA_noXY_samp[,1])
unique(noSex[,1])
table(miRNA_noXY_samp[,1])
# chr1  chr10  chr11  chr12  chr13  chr14  chr15  chr16  chr17  chr18  chr19   chr2  chr20  chr21  chr22 
#  81     26     45     35     19     23     22     27     62      5     55     34     12      9     34 
# chr3   chr4   chr5   chr6   chr7   chr8   chr9 chrEBV 
#  38     17     30     37     40     26     39      6 

length(grep(pattern = "novel", x = rownames(miRNA_noXY_samp), value = TRUE))
#149  novel miRNAs persist
idx_novel <- grep(pattern = "novel", x = rownames(miRNA_noXY_samp), value = FALSE)
length(idx_novel)

Novel_miRNA_pass_exp <- rownames(miRNA_noXY_samp)[idx_novel]
length(Novel_miRNA_pass_exp) # 149
saveRDS(Novel_miRNA_pass_exp, "Novel_miRNA_pass_exp.rds")

summary(miRNA_noXY_samp[,3])

# plotting the novel miRNA expression levels
melt_novel <- melt(miRNA_noXY_samp[idx_novel,-1])
head(melt_novel)
ggplot(data = melt_novel, mapping = aes(y=value, x=variable))+
  geom_boxplot()

library(reshape)
head(idx_novel)
head(miRNA_noXY_samp)
dat<-log(miRNA_noXY_samp[idx_novel,-1]+1)
#  novel_mean<-apply(dat, MARGIN = 1, mean)
#  boxplot(novel_mean, ylim=c(0,500))
length(idx_novel)
dat$group <- row.names(dat)
nrow(dat) #74 pass exp thres
dat.m <- melt(dat, id.vars = "group")
# made pdf Expression_of_novel_miRNAs.pdf
ggplot(dat.m, aes(group, value)) + geom_boxplot() +
  ylab("log(count+1)") +
  xlab("Novel miRNA")

# filter out EBV miRNAs because of Copy Number skewing
miRNA_noXY_samp_noEBV <- miRNA_noXY_samp[miRNA_noXY_samp$Chr != "chrEBV",]
head(miRNA_noXY_samp_noEBV)
dim(miRNA_noXY_samp_noEBV) # 716  18

filter_miRNA_counts <- readRDS("../filter_miRNA_counts.rds")
idx_novel_fil <- which(rownames(miRNA_noXY_samp_noEBV) %in% filter_miRNA_counts)
length(idx_novel_fil) #9
miRNA_noXY_samp_noEBV_fil <- miRNA_noXY_samp_noEBV[-idx_novel_fil,]
dim(miRNA_noXY_samp_noEBV_fil) #707  18
saveRDS(miRNA_noXY_samp_noEBV_fil, "miRNA_noXYEBV_filt.rds")
head(miRNA_noXY_samp_noEBV_fil)

# make a bed file for the final miRNA count set so we can intersect with haplotypes
Final_miRNA_count_set <- rownames(miRNA_noXY_samp_noEBV_fil)
head(Final_miRNA_count_set)
saveRDS(Final_miRNA_count_set, "Final_miRNA_count_set.rds")
miRNA_count_file <- read.table("GM77_miRNA.txt", header=F,
                        sep="\t")
str(miRNA_count_file)
idx_miRNA_final_count <- which(miRNA_count_file$V1 %in% Final_miRNA_count_set)
length(idx_miRNA_final_count) # 707
miRNA_count_file_final <- miRNA_count_file[idx_miRNA_final_count,]
head(miRNA_count_file_final)
dim(miRNA_count_file_final) #707
miRNA_count_final_bed <- miRNA_count_file_final[,c(2,3,4,1,6,5)]
write.table(miRNA_count_final_bed, file="miRNA_count_final.bed", append = F, 
            quote = F, sep = "\t", row.names = F, col.names = F)

## TO-Do, PCA the raw miRNAs
# Do EDA seq of them, and look at UQ norm/or TMM norm
# Do PCA heatmaps, correlation, etc.
# 
head(miRNA_noXY_samp_noEBV_fil)
colnames(miRNA_noXY_samp_noEBV_fil)[-1] <- substr(colnames(miRNA_noXY_samp_noEBV_fil)[-1], start=4, stop = 7)
facts <- as.factor(colnames(miRNA_noXY_samp_noEBV_fil)[-1])
facts
set <- newSeqExpressionSet(as.matrix(miRNA_noXY_samp_noEBV_fil[,-1]),
                           phenoData = data.frame(facts, 
                                                  row.names=colnames(miRNA_noXY_samp_noEBV_fil)[-1]))

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11, "Spectral"))(17)
color_sex <- brewer.pal(3, "Set1")
sex <- batch$sex

pdf(file = "RLE_miRNA_noXY_noEBV_fil.pdf", width = 6, height = 4)
plotRLE(set, outline=FALSE, ylim=c(-3, 3), col=colors[facts]) 
title("miRNA b4 Normalization")
dev.off()

pdf(file = "PCA_miRNA_noXY_noEBV_fil.pdf", width = 6, height = 6)
plotPCA(set, col=colors[facts], cex=1.2, k=2)
plotPCA(set, col=colors[facts], cex=1.2, k=3)
plotPCA(set, col=color_sex[sex], cex=1.2, k=2)
plotPCA(set, col=color_sex[sex], cex=1.2, k=3)
dev.off()

# perform TMM normalization
set_DGE <- DGEList(counts = set@assayData$counts, group = facts)
set_norm <- calcNormFactors(set_DGE, method = "TMM")
# it's fine to use CPM for between sample comparison
count_norm <- cpm(set_norm, normalized.lib.sizes=TRUE)
head(count_norm)
saveRDS(count_norm, "miRNA_noXYEBV_filt_TMM.rds")
# load the TMM nomalized gene counts back into EDAseq
set_TMM <- newSeqExpressionSet(as.matrix(count_norm),
                               phenoData = data.frame(facts, 
                                                      row.names=colnames(count_norm)))

pdf(file = "RLE_miRNA_noXY_noEBV_fil_TMM.pdf", width = 6, height = 4)
plotRLE(set_TMM, outline=FALSE, ylim=c(-3, 3), col=colors[facts]) 
title("miRNA noSexEbv filt TMM")
dev.off()

pdf(file = "PCA_miRNA_noXY_noEBV_fil_TMM.pdf", width = 6, height = 6)
plotPCA(set_TMM, col=colors[facts], cex=1.2, k=2)
title("miRNA noXYEBV filt TMM")
plotPCA(set_TMM, col=colors[facts], cex=1.2, k=3)
plotPCA(set_TMM, col=color_sex[sex], cex=1.2, k=2)
title("miRNA noXYEBV filt TMM")
plotPCA(set_TMM, col=color_sex[sex], cex=1.2, k=3)
dev.off()

# I checked upper quantile normalization but it wasn't as effective as TMM
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[facts]) 
title("RLE plot after UQ norm")
plotPCA(set, col=colors[facts], cex=1.2, k=3)

head(count_norm)
dim(count_norm)

# Normalize by total Reads
## PCA Analysis

counts_load <- readRDS("miRNA_noXYEBV_filt_TMM.rds")
head(counts_load)
dim(counts_load)

x_raw_noXYEBV_fil <- princomp(miRNA_noXY_samp_noEBV_fil[,-1])
x_raw_noXYEBV_fil_TMM <- princomp(counts_load)
####PCA Pvalue plot 
batch <- read.table("../../metadata.txt", header=T)
head(batch)
str(batch)

batch[,7] <- as.factor(batch[,7])
batch[,8] <- as.factor(batch[,8])
batch[,9] <- as.factor(batch[,9])
colnames(batch)[7] <- "lane"
colnames(batch)[9] <- "gen."

# add EBV
EBV_covariate_value<- readRDS("../../../../../../RNA-seq/Sleuth/EBV_covariate_value.rds")

batch$EBV <- cut(EBV_covariate_value, breaks=4)

#pca <- x_raw
pca <- x_raw_noXYEBV_fil
pca <- x_raw_noXYEBV_fil_TMM

Sex.p<-c()
for (i in seq(1:10)){
  Sex<-batch[,3]
  Sex<-lm(pca$loadings[,i]~Sex)
  Sex.p<-c(Sex.p, anova(Sex)$Pr[1])
}

Mito.p<-c()
for (i in seq(1:10)){
  Mito<-batch[,4]
  Mito<-lm(pca$loadings[,i]~Mito)
  Mito.p<-c(Mito.p, anova(Mito)$Pr[1])
}

harvest.p<-c()
for (i in seq(1:10)){
  harvest<-batch[,5]
  harvest<-lm(pca$loadings[,i]~harvest)
  harvest.p<-c(harvest.p, anova(harvest)$Pr[1])
}

extract.p<-c()
for (i in seq(1:10)){
  extract<-batch[,6]
  extract<-lm(pca$loadings[,i]~extract)
  extract.p<-c(extract.p, anova(extract)$Pr[1])
}

Lane1.p<-c()
for (i in seq(1:10)){
  Lane1<-batch[,7]
  Lane1<-lm(pca$loadings[,i]~Lane1)
  Lane1.p<-c(Lane1.p, anova(Lane1)$Pr[1])
}

generation.p<-c()
for (i in seq(1:10)){
  generation<-batch[,9]
  generation<-lm(pca$loadings[,i]~generation)
  generation.p<-c(generation.p, anova(generation)$Pr[1])
}

EBV.p<-c()
for (i in seq(1:10)){
  EBV<-batch[,10]
  EBV<-lm(pca$loadings[,i]~EBV)
  EBV.p<-c(EBV.p, anova(EBV)$Pr[1])
}

pvals.raw <- rbind(Sex.p,Mito.p,harvest.p,extract.p,Lane1.p, generation.p,EBV.p)
pvals.raw <- data.matrix(pvals.raw)
rownames(pvals.raw) <- colnames(batch[,-c(1,2,8)])
library(RColorBrewer)
library(gplots)
hmcol4<-colorRampPalette(brewer.pal(9,"Blues"))(5)
hmcol4[1]<-"white"
logpvals.raw <- -log10(pvals.raw)
max <- ifelse(max(logpvals.raw)>5, max(logpvals.raw), 5)
# pdf(file = "PCA_heatmap_miRNA_noXY_noEBV_fil.pdf", width = 9, height = 6)
pdf(file = "PCA_heatmap_miRNA_noXY_noEBV_fil_TMM2.pdf", width = 9, height = 6)
heatmap.2(logpvals.raw,Rowv=F,Colv=colnames(pvals.raw),dendrogram='none',trace='none',
          margins=c(8,8),colsep=c(1:11),rowsep=c(1:16),
          sepwidth=c(0.025,0.025), sepcolor="black", col=hmcol4, 
          breaks=c(0,1.30103,2,3,4,max),
          key.xlab = "-log10 pvalue", main="miRNA noSexEBV filt TMM")  
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
# pdf(file = "PCA_varBar_miRNA_noXY_noEBV_fil.pdf", width = 6, height = 4)
pdf(file = "PCA_varBar_miRNA_noXY_noEBV_fil_TMM.pdf", width = 6, height = 4)
barplot(100*p.variance.explained10, las=2, xlab='', ylab='% Variance Explained', 
        main="Variance Explained by component", 
        ylim = c(0,100))
dev.off()
