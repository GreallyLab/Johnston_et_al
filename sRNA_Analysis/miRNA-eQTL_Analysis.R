setwd("/Volumes/users-1/Andrew Johnston/17_member/QTLs_2017/miRNA_QTLs/")
# loading necessary libraries
library(qvalue)
library(ggplot2)
library(gridExtra)
library(Rmisc)
options(stringsAsFactors=F)

# making an object with the sex of the children
Sex_child <- c("F","F","F","M","M","M","F","M","F","M","M")
length(Sex_child)
Sex_child_factor <- as.factor(Sex_child)

# Loading in the Haplotype data
Haplo_data <- read.table("../haplotype_transmission_grch38_named_noXYHeader.txt")
head(Haplo_data)
dim(Haplo_data) # [1] 708  17
# changing the format from AC, AD, BC, or BD to paternal 0 or 1 and maternal 0 or 1
for (i in 1:nrow(Haplo_data)){
  temp_haplo <- rep(NA, 22)
  k<-1
  for (j in 6:16) {
    if (Haplo_data[i,j] == "AC"){
      temp_haplo[k] <- 0
      k <- k+1
      temp_haplo[k] <- 0
      k <- k+1
    }
    if (Haplo_data[i,j] == "AD"){
      temp_haplo[k] <- 0
      k <- k+1
      temp_haplo[k] <- 1
      k <- k+1
    }
    if (Haplo_data[i,j] == "BC"){
      temp_haplo[k] <- 1
      k <- k+1
      temp_haplo[k] <- 0
      k <- k+1
    }
    if (Haplo_data[i,j] == "BD"){
      temp_haplo[k] <- 1
      k <- k+1
      temp_haplo[k] <- 1
      k <- k+1
    }
  }
  if (i == 1){
    temp_haplo_df <- temp_haplo
  }
  if (i > 1){
    temp_haplo_df <- rbind(temp_haplo_df,temp_haplo)
  }
}
head(temp_haplo_df)
rownames(temp_haplo_df) <- rownames(Haplo_data)
Haplo_data_edit<- cbind(Haplo_data,temp_haplo_df)

# removing the columns that I no longer need
Haplo_data_edit <- Haplo_data_edit[,-(4:16)]
head(Haplo_data_edit)

# loading in the transcripts intersected with platinum haplotype information
miRNA_haplo_inter <- read.table("../miRNA_haplo_intersect.txt")
head(miRNA_haplo_inter)
table(miRNA_haplo_inter$V1)  # no XY chr
table(table(miRNA_haplo_inter$V4)) # 707 miRNA tested;
dim(miRNA_haplo_inter) # 707   11

# Loading in the normalized counts of transcripts as the quantitative trait
# TMM norm
miRNA_vals_TMM <- as.data.frame(readRDS("../../miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Novel_miRNA/Count_miRNA_novel/miRNA_noXYEBV_filt_TMM.rds"))
dim(miRNA_vals_TMM) #707    17
head(miRNA_vals_TMM)
miRNA_vals_TMM$Gene_name <- row.names(miRNA_vals_TMM)

# merging the miRNA count values with the haplotype peak intersection information
miRNA_inter_vals <- merge(x=miRNA_haplo_inter, y=miRNA_vals_TMM, by.x = "V4", by.y = "Gene_name", all.y = TRUE,sort = F)
head(miRNA_inter_vals)
dim(miRNA_inter_vals) # 707 28 

head(Haplo_data_edit)
# merging the miRNA count values / haplotype-peak intersect with haplo-genotype information
miRNA_haplo <- merge(x = miRNA_inter_vals, y=Haplo_data_edit, by.x = "V10", by.y = "V17", all.x=T, sort = F)
head(miRNA_haplo)
dim(miRNA_haplo) # [1] 707     53

#setup the objects for anova loop
miRNA_names <- miRNA_haplo$V4
length(miRNA_names) # 707
counts<- miRNA_haplo[,c(14:23,28)]
head(counts)
dim(counts)
miRNA_pat <- miRNA_haplo[,c(32,34,36,38,40,42,44,46,48,50,52)]
head(miRNA_pat)
dim(miRNA_pat)
miRNA_mat <- miRNA_haplo[,c(33,35,37,39,41,43,45,47,49,51,53)]
head(miRNA_mat)
dim(miRNA_mat)

# Run a loop which computes an anova for each test: 
# Is the model of paternal/maternal haplotype affecting miRNA counts better than the null?
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

i<- which(names(miRNA_anova_pval) == "MIMAT0027384")
plot_count_haplo_list[[6]]

miRNA_anova_pval <- rep(1,nrow(counts))
for (i in 1:nrow(counts)){
  y <- t(counts[i,])
  m <- t(miRNA_mat[i,])
  p <- t(miRNA_pat[i,])
  temp_df <- cbind(y,m,p)
  colnames(temp_df) <- c("y","m","p")
  temp_df
  class(temp_df)
  full_model <- lm(formula = y ~ m + p, data=as.data.frame(temp_df))
  if (anova(full_model)$Pr[1] < anova(full_model)$Pr[2]){
    anova_p <- anova(full_model)$Pr[1]
  }
  else {
    anova_p <- anova(full_model)$Pr[2]
  }
  miRNA_anova_pval[i] <- anova_p
}
head(miRNA_anova_pval)
length(miRNA_anova_pval) # 707
sum(is.na(miRNA_anova_pval)) # 0
names(miRNA_anova_pval)<-miRNA_names
sum(miRNA_anova_pval==0) # 0
which.min(miRNA_anova_pval)
miRNA_anova_pval[340] #novelMiR_268 : 8.219469e-07

# plot distribution of pvalues 
pdf(file="qplot_miRNA-eQTL_p-value.pdf", width=6, height = 4)
qplot(miRNA_anova_pval, main="miRNA-eQTL p-value distribution", xlab="p-value") 
dev.off()
summary(miRNA_anova_pval)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.1317000 0.2838000 0.3273000 0.4861000 0.9353000 

# run qvalue (BJH) correction to the pvalues to get "significant" miRNA-eQTLs
library(qvalue)
qvals_miRNA_anova <- qvalue(miRNA_anova_pval, fdr.level=.05)
head(qvals_miRNA_anova$qvalues)
sum(qvals_miRNA_anova$significant, na.rm = T) # 19 vs. 37-TMM
sum(qvals_miRNA_anova$qvalues < .4,  na.rm = T) # 18 vs. 12-TMM
sum(qvals_miRNA_anova$qvalues < .05,  na.rm = T) # 2 vs. 1
sum(qvals_miRNA_anova$qvalues < .017,  na.rm = T) # 2 vs. 1


# looking at distribution of the qvalues
summary(qvals_miRNA_anova$qvalues)
#       Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.001777 0.773000 0.814100    0.782800  0.839600  0.873700 

# plotting the distribution of the chrQTLs and placing a red line at the cutoff
pdf(file="qplot_miRNA-eQTL_q-value.pdf", width=6, height = 4)
qplot(qvals_miRNA_anova$qvalues, main="miRNA-eQTL q-value distribution TMM", xlab="q-value") +
  geom_vline(xintercept = .5, colour="red", linetype = "longdash")
dev.off()

qvals_miRNA_anova_TMM <- qvals_miRNA_anova
saveRDS(qvals_miRNA_anova_TMM, "qvals_miRNA_anova_TMM.rds")
qvals_miRNA_anova_TMM <- readRDS("qvals_miRNA_anova_TMM.rds")
qvals_miRNA_anova <- qvals_miRNA_anova_TMM
head(qvals_miRNA_anova$qvalues)

# using .5 as the qvalue cutoff for the miRNA-eQTLs
head(qvals_miRNA_anova$qvalues)
idx_sig_miRNA_names <- names(qvals_miRNA_anova$qvalues)[qvals_miRNA_anova$significant]
length(idx_sig_miRNA_names) #37

# get the index for the significant peaks
idx_sig_miRNA <- which(qvals_miRNA_anova$qvalues < .5 & !is.na(qvals_miRNA_anova$qvalues))
length(idx_sig_miRNA) #37

idx_sig_miRNA <- which(names(qvals_miRNA_anova$qvalues) %in% cat_p_m)
length(idx_sig_miRNA)

idx_sig_miRNA_names <- names(qvals_miRNA_anova$qvalues)[idx_sig_miRNA]

# plot the signficant peaks - .5
plot_count_haplo_list <- list()
i<-13
f<-1
for (i in idx_sig_miRNA){
  y <- t(counts[i,])
  m <- t(miRNA_mat[i,])
  p <- t(miRNA_pat[i,])
  temp_df <- cbind(y,m,p, 0,Sex_child_factor)
  colnames(temp_df) <- c("y","m","p","g","sex")
  temp_df <- as.data.frame(temp_df)
  for (j in 1:nrow(temp_df)) {
    if (temp_df$m[j] == 0 && temp_df$p[j] == 0){temp_df$g[j] <- "0/0"}
    if (temp_df$m[j] == 1 && temp_df$p[j] == 0){temp_df$g[j] <- "1/0"}
    if (temp_df$m[j] == 0 && temp_df$p[j] == 1){temp_df$g[j] <- "0/1"}
    if (temp_df$m[j] == 1 && temp_df$p[j] == 1){temp_df$g[j] <- "1/1"}
  }
  temp_plot <- ggplot(temp_df, aes(factor(g), y))
  temp_plot <- temp_plot + geom_boxplot() + 
    geom_jitter(aes(colour=factor(sex,labels=c("Female","Male")))) +
    scale_color_discrete(name = "Gender") +
    ggtitle(miRNA_haplo$V4[i]) +
    ylab("miRNA Counts") + xlab("Haplotype")
  plot_count_haplo_list[[f]] <- temp_plot
  f <- f+1
}
length(plot_count_haplo_list) # 47
random_samp <- sample(seq(1,19), 6)
multiplot(plotlist = plot_count_haplo_list[random_samp], cols = 2)

# chosen miRNA QTLs for poster (Marmur)
multiplot(plotlist = plot_count_haplo_list[c(15,22)], cols = 2)


plot_count_haplo_list[[20]]
tail(seq(6,42,6))
pdf(file = "miRNA_genotype_eQTL_.5.pdf")
j<-1
tail(seq(6,42,6))
for (i in seq(6,42,6)){
  multiplot(plotlist = plot_count_haplo_list[j:i], cols = 2)
  j <- i+1
}
dev.off()
saveRDS(plot_count_haplo_list, "plot_count_haplo_list.rds")
rm(plot_count_haplo_list)
plot_count_haplo_list2 <- readRDS("plot_count_haplo_list.rds")
plot_count_haplo_list2[[5]]


# generate a table of information for the QTL miRNAs
miRNA_haplo$eQTL_qval <- qvals_miRNA_anova$qvalues
miRNA_haplo_eQTL <- miRNA_haplo[idx_sig_miRNA,]

dim(miRNA_haplo_eQTL)
head(miRNA_haplo_eQTL)
colnames(miRNA_haplo_eQTL)
miRNA_haplo_eQTL_edit <- miRNA_haplo_eQTL[,c(3,4,5,2,7,6,1,(12:28),54)]
head(miRNA_haplo_eQTL_edit)
saveRDS(miRNA_haplo_eQTL_edit, "miRNA_haplo_eQTL_edit.rds")
miRNA_haplo_eQTL_edit <- readRDS("miRNA_haplo_eQTL_edit.rds")
head(miRNA_haplo_eQTL_edit)
sum(idx_sig_miRNA_names %in% miRNA_haplo_eQTL_edit$V4) #13
miRNA_haplo_eQTL_edit[!(miRNA_haplo_eQTL_edit$V4 %in% idx_sig_miRNA_names), 4] #13
miRNA_haplo_eQTL_edit$V4
idx_sig_miRNA_names[!(idx_sig_miRNA_names %in% miRNA_haplo_eQTL_edit$V4)]
