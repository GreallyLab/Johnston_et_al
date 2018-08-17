setwd("/Volumes/users-1/Andrew Johnston/17_member/QTLs_2017/ATAC_summit_QTL_2/")
# loading necessary libraries
library(qvalue)
library(ggplot2)
library(gridExtra)
library(Rmisc)

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
saveRDS(Haplo_data_edit, "../Haplo_data_edit_forQTL.rds")
sum(is.na(Haplo_data_edit)) #0

# loading in the ATACseq summits intersected with platinum haplotype information
ATAC_haplo_inter <- read.table("ATAC_summit_Plat_haplo_interesect.txt")
head(ATAC_haplo_inter)
table(ATAC_haplo_inter$V1)  # no sex chr
table(table(ATAC_haplo_inter$V4)) # 139301 summits tested
dim(ATAC_haplo_inter) # 139301     11
sum(is.na(ATAC_haplo_inter)) #0

# Loading in the mean values of ATACseq peaks as the quantitative trait
ATAC_vals <- read.table("../../ATACseq/Quant_summit_2/Summit_counts_CQN_norm_child.txt", header=T)
dim(ATAC_vals) #57407     35
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
                             gm93=(ATAC_vals[,33]+ATAC_vals[,34])/2,
                             Peak_number=ATAC_vals[,35])
head(ATAC_vals_mean)
head(ATAC_haplo_inter)


# merging the ATAC peak values with the haplotype peak intersection information
ATAC_inter_vals <- merge(x=ATAC_haplo_inter, y=ATAC_vals_mean, by.x = "V4", by.y = "Peak_number",sort = F)
head(ATAC_inter_vals)
dim(ATAC_inter_vals) # 57314     28

# merging the ATAC peak values / haplotype-peak intersect with haplo-genotype information
ATAC_haplo <- merge(x = ATAC_inter_vals, y=Haplo_data_edit, by.x = "V10", by.y = "V17", all.x=T, sort = F)
head(ATAC_haplo)
dim(ATAC_haplo) # [1] 57314     53
sum(is.na(ATAC_haplo)) #0

#setup the objects for anova loop
peaks <- ATAC_haplo$V4
length(peaks)
counts<- ATAC_haplo[,c(14:23,28)]
head(counts)
dim(counts)
sum(is.na(counts))
ATAC_pat <- ATAC_haplo[,c(32,34,36,38,40,42,44,46,48,50,52)]
head(ATAC_pat)
dim(ATAC_pat)
sum(is.na(ATAC_pat)) #0
ATAC_mat <- ATAC_haplo[,c(33,35,37,39,41,43,45,47,49,51,53)]
head(ATAC_mat)
dim(ATAC_mat)
sum(is.na(ATAC_mat)) #0

# Run a loop which computes an anova for each test: 
# Is the model of paternal/maternal haplotype affecting ATAC peak strength better than the null?
ATAC_anova_pval_m <- rep(1,nrow(counts))
ATAC_anova_pval_p <- rep(1,nrow(counts))

for (i in 1:nrow(counts)){
  y <- t(counts[i,])
  m <- t(ATAC_mat[i,])
  p <- t(ATAC_pat[i,])
  temp_df <- cbind(y,m,p)
  colnames(temp_df) <- c("y","m","p")
  full_model <- lm(formula = y ~ m + p, data=as.data.frame(temp_df))
  ATAC_anova_pval_m[i] <- anova(full_model)$Pr[1] 
  ATAC_anova_pval_p[i] <- anova(full_model)$Pr[2]
}
# In anova.lm(full_model) :
# ANOVA F-tests on an essentially perfect fit are unreliable
# I got 6 warnings of this nature.. are there perfect fits?

length(ATAC_anova_pval_m) # 57314
sum(is.na(ATAC_anova_pval_m)) # 0
names(ATAC_anova_pval_m)<-peaks
sum(ATAC_anova_pval_m==0) # 0
pdf(file="qplot_chrQTL_p-value_mat.pdf", width=6, height = 4)
qplot(ATAC_anova_pval_m, main="chrQTL maternal p-value distribution", xlab="p-value") 
dev.off()
summary(ATAC_anova_pval_m)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000006 0.2572000 0.5073000 0.5051000 0.7561000 1.0000000 
library(qvalue)
qvals_ATAC_anova_m <- qvalue(ATAC_anova_pval_m, fdr.level=.5)
head(qvals_ATAC_anova_m$qvalues)
sum(qvals_ATAC_anova_m$significant, na.rm = T) # 74
sum(qvals_ATAC_anova_m$qvalues < .4,  na.rm = T) # 41
sum(qvals_ATAC_anova_m$qvalues < .05,  na.rm = T) # 5
pdf(file="qplot_chrQTL_q-value_mat.pdf", width=6, height = 4)
qplot(qvals_ATAC_anova_m$qvalues, main="chrQTL maternal q-value distribution", xlab="q-value") +
  geom_vline(xintercept = .5, colour="red", linetype = "longdash")
dev.off()

length(ATAC_anova_pval_p) # 57314
sum(is.na(ATAC_anova_pval_p)) # 43
names(ATAC_anova_pval_p)<-peaks
ATAC_anova_pval_p_noNA <- ATAC_anova_pval_p[!is.na(ATAC_anova_pval_p)]
length(ATAC_anova_pval_p_noNA) #57271
sum(ATAC_anova_pval_p_noNA==0) # 0
pdf(file="qplot_chrQTL_p-value_pat.pdf", width=6, height = 4)
qplot(ATAC_anova_pval_p_noNA, main="chrQTL paternal p-value distribution", xlab="p-value") 
dev.off()
summary(ATAC_anova_pval_p_noNA)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000001 0.2198000 0.4697000 0.4782000 0.7303000 1.0000000 
library(qvalue)
qvals_ATAC_anova_p <- qvalue(ATAC_anova_pval_p_noNA, fdr.level=.5)
head(qvals_ATAC_anova_p$qvalues)
sum(qvals_ATAC_anova_p$significant, na.rm = T) #  216
sum(qvals_ATAC_anova_p$qvalues < .4,  na.rm = T) # 129
sum(qvals_ATAC_anova_p$qvalues < .05,  na.rm = T) # 12
pdf(file="qplot_chrQTL_q-value_pat.pdf", width=6, height = 4)
qplot(qvals_ATAC_anova_p$qvalues, main="chrQTL paternal q-value distribution", xlab="q-value") +
  geom_vline(xintercept = .5, colour="red", linetype = "longdash")
dev.off()

idx_sig_ATAC_names_m <- names(qvals_ATAC_anova_m$qvalues)[qvals_ATAC_anova_m$significant]
idx_sig_ATAC_names_p <- names(qvals_ATAC_anova_p$qvalues)[qvals_ATAC_anova_p$significant]
#length(idx_sig_ATAC_names) #19
sum(idx_sig_ATAC_names_m %in% idx_sig_ATAC_names_p) #14 (so it should be 276?)
idx_sig_ATAC_names_m[idx_sig_ATAC_names_m %in% idx_sig_ATAC_names_p] 

cat_p_m <- c(idx_sig_ATAC_names_p, idx_sig_ATAC_names_m)
sum(duplicated(cat_p_m)) # 14

# get the index for the significant peaks
idx_sig_peaks <- which(peaks %in% cat_p_m)
length(idx_sig_peaks) # 276

idx_sig_ATAC_names <- names(qvals_ATAC_anova_m$qvalues)[idx_sig_ATAC]


# plot the signficant peaks - .5
plot_count_haplo_list <- list()
f <- 1
i <-457
for (i in idx_sig_peaks){
  y <- t(counts[i,])
  m <- t(ATAC_mat[i,])
  p <- t(ATAC_pat[i,])
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
    ggtitle(ATAC_haplo$V4[i]) +
    ylab("Peak_strength") + xlab("Haplotype")
  plot_count_haplo_list[[f]] <- temp_plot
  f <- f+1
}
length(plot_count_haplo_list) # 276
random_samp <- sample(seq(1,276), 6)
multiplot(plotlist = plot_count_haplo_list[random_samp], cols = 2)
plot_count_haplo_list[[240]]

tail(seq(6,276,6))
pdf(file = "ATAC_summit_genotype_QTL_.5.pdf")
j<-1
for (i in seq(6,276,6)){
  multiplot(plotlist = plot_count_haplo_list[j:i], cols = 2)
  j <- i+1
}
dev.off()
#multiplot(plotlist = plot_count_haplo_list[1:6], cols = 2)
#multiplot(plotlist = plot_count_haplo_list[6:12], cols = 2)
#multiplot(plotlist = plot_count_haplo_list[c(12,22)], cols = 2)

saveRDS(plot_count_haplo_list,"plot_count_haplo_list_child_summit.rds")
plot_count_haplo_list <- readRDS("plot_count_haplo_list_child_summit.rds")
rm(plot_count_haplo_list)

idx_same_matpat <- which(idx_sig_ATAC_names_m %in% idx_sig_ATAC_names_p)

#matpat_sig_qvalue <- qvals_ATAC_anova_p
pat_sig_qvalue <- qvals_ATAC_anova_p$qvalues[names(qvals_ATAC_anova_p$qvalues) %in%
                                                idx_sig_ATAC_names_p]
mat_sig_qvalue <- qvals_ATAC_anova_m$qvalues[names(qvals_ATAC_anova_m$qvalues) %in%
                                                idx_sig_ATAC_names_m]
pat_sig_qvalue_df <- data.frame(name=names(pat_sig_qvalue), qvalue=pat_sig_qvalue)
mat_sig_qvalue_df <- data.frame(name=names(mat_sig_qvalue), qvalue=mat_sig_qvalue)
cat_matpat_df <- rbind(mat_sig_qvalue_df, pat_sig_qvalue_df)
rownames(cat_matpat_df) <- seq(1,nrow(cat_matpat_df),1)
dup_matpat <- cat_matpat_df$name[duplicated(cat_matpat_df$name)]
i<-1
for (i in 1:length(dup_matpat)) {
  dups <- cat_matpat_df[cat_matpat_df$name == dup_matpat[i],]
  if(dups[1,2]<dups[2,2]) {
    cat_matpat_df <- cat_matpat_df[-as.numeric(rownames(dups)[2]),]
  }
  else {
    cat_matpat_df <- cat_matpat_df[-as.numeric(rownames(dups)[1]),]
  }
  rownames(cat_matpat_df) <- seq(1,nrow(cat_matpat_df),1)
}
cat_matpat_df

ATAC_haplo_chrQTL <- merge(ATAC_haplo, y=cat_matpat_df, by.x="V4", by.y="name")
dim(ATAC_haplo_chrQTL) #276 54
head(ATAC_haplo_chrQTL)

colnames(ATAC_haplo_chrQTL)
ATAC_haplo_chrQTL_edit <- ATAC_haplo_chrQTL[,c(3,4,5,1,6,7,2,(12:28),54)]
head(ATAC_haplo_chrQTL_edit)
saveRDS(ATAC_haplo_chrQTL_edit, "ATAC_haplo_chrQTL_edit_matpat.rds")

which.min(ATAC_haplo_chrQTL_edit$qvalue)
ATAC_haplo_chrQTL_edit[59,]
which(idx_sig_ATAC_names %in% ATAC_haplo_chrQTL_edit[59,4])

### getting the mat, pat , or both for each summit 
head(cat_p_m)

idx_sig_ATAC_names_m <- names(qvals_ATAC_anova_m$qvalues)[qvals_ATAC_anova_m$significant]
head(idx_sig_ATAC_names_m)
df_sig_ATAC_names_m <- data.frame(Parent="Mat", Summit=idx_sig_ATAC_names_m)
head(df_sig_ATAC_names_m)
idx_sig_ATAC_names_p <- names(qvals_ATAC_anova_p$qvalues)[qvals_ATAC_anova_p$significant]
df_sig_ATAC_names_p <- data.frame(Parent="Pat", Summit=idx_sig_ATAC_names_p)
head(df_sig_ATAC_names_p)

df_df_sig_ATAC_names <- merge(df_sig_ATAC_names_m, df_sig_ATAC_names_p, by="Summit", all=TRUE)
head(df_df_sig_ATAC_names) 

for (i in 1:nrow(df_df_sig_ATAC_names)){
  if (!is.na(df_df_sig_ATAC_names$Parent.x[i]) &
      !is.na(df_df_sig_ATAC_names$Parent.y[i])) {
    df_df_sig_ATAC_names$Allele[i] <- "Both"
  }
  if (is.na(df_df_sig_ATAC_names$Parent.x[i])){
    df_df_sig_ATAC_names$Allele[i] <- "Pat"
  }
  if (is.na(df_df_sig_ATAC_names$Parent.y[i])){
    df_df_sig_ATAC_names$Allele[i] <- "Mat"
  }
}
head(df_df_sig_ATAC_names)
table(df_df_sig_ATAC_names$Allele)
df_sig_ATAC_names <- df_df_sig_ATAC_names[,c(1,4)]
head(df_sig_ATAC_names)
saveRDS(df_sig_ATAC_names,"df_sig_ATAC_names.rds")

#length(idx_sig_ATAC_names) #19
sum(idx_sig_ATAC_names_m %in% idx_sig_ATAC_names_p) #14 (so it should be 276?)
idx_sig_ATAC_names_m[idx_sig_ATAC_names_m %in% idx_sig_ATAC_names_p] 

cat_p_m <- c(idx_sig_ATAC_names_p, idx_sig_ATAC_names_m)

# Getting the values for the chrQTLs for each summit
chrQTL_value_list <-list()
f<-1
for (i in idx_sig_peaks){
  y <- t(counts[i,])
  m <- t(ATAC_mat[i,])
  p <- t(ATAC_pat[i,])
  temp_df <- cbind(y,m,p, 0,Sex_child_factor)
  colnames(temp_df) <- c("y","m","p","g","sex")
  temp_df <- as.data.frame(temp_df)
  for (j in 1:nrow(temp_df)) {
    if (temp_df$m[j] == 0 && temp_df$p[j] == 0){temp_df$g[j] <- "0/0"}
    if (temp_df$m[j] == 1 && temp_df$p[j] == 0){temp_df$g[j] <- "1/0"}
    if (temp_df$m[j] == 0 && temp_df$p[j] == 1){temp_df$g[j] <- "0/1"}
    if (temp_df$m[j] == 1 && temp_df$p[j] == 1){temp_df$g[j] <- "1/1"}
  }
  temp_df <- temp_df[,c(1,4,5)]
  chrQTL_value_list[[f]] <- temp_df
  names(chrQTL_value_list)[f] <- as.character(ATAC_haplo$V4[i])
  f<-f+1
}
length(chrQTL_value_list)

saveRDS(chrQTL_value_list,"chrQTL_value_list.rds")
