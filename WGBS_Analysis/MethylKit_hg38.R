#' this is a test run for the methylKit script
#+ load-libs, include = TRUE
suppressPackageStartupMessages({
	library(methylKit)
	library(RColorBrewer)
	library(gplots)
	library(knitr)
})

#+ setup, include=FALSE
opts_chunk$set(fig.width = 7, fig.height = 7, fig.path='Figs_hg38_2/')

#+ setup-functions, include=FALSE
generate_heatmap <- function(x,y){
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  cormat <- round(cor(x),3)
  melted_cormat <- melt(cormat)
  # before reordering
  temp_gg <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + ggtitle(label=paste(y,"correlation",sep = " ")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))
  print(temp_gg)
  cormat_noclust <- cormat
  cormat <- reorder_cormat(cormat)
  melted_cormat <- melt(cormat, na.rm = TRUE)
  temp_gg <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + ggtitle(label=paste(y,"clustered correlation",sep = " ")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))
  print(temp_gg)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()+
    ggtitle(label=paste(y,"Clustered Correlation Plot",sep = " "))
  
  # Print the heatmap
  print(ggheatmap)
  cormat_noclust
}
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

#' Create file list
#+ create-file-list, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60), results="hide"
file.list=list("gm77/analysis_hg38/BS-TAG-Sample-77_CpG_percentFix.methylKit",
               "gm78/analysis_hg38/BS-TAG-Sample-78_CpG_percentFix.methylKit",
               "gm79/analysis_hg38/BS-TAG-Sample-79_CpG_percentFix.methylKit",
               "gm80/analysis_hg38/BS-TAG-Sample-80_CpG_percentFix.methylKit",
               "gm81/analysis_hg38/BS-TAG-Sample-81_CpG_percentFix.methylKit",
               "gm82/analysis_hg38/BS-TAG-Sample-82_CpG_percentFix.methylKit",
               "gm83/analysis_hg38/BS-TAG-Sample-83_CpG_percentFix.methylKit",
               "gm84/analysis_hg38/BSTag-Sample-84-reprep_CpG_percentFix.methylKit",
               "gm85/analysis_hg38/BS-TAG-Sample-85_CpG_percentFix.methylKit",
               "gm86/analysis_hg38/BS-TAG-Sample-86_CpG_percentFix.methylKit",
               "gm87/analysis_hg38/BS-TAG-Sample-87_CpG_percentFix.methylKit",
               "gm88/analysis_hg38/BS-TAG-Sample-88_CpG_percentFix.methylKit",
               "gm89/analysis_hg38/BSTag-Sample-89-reprep_CpG_percentFix.methylKit",
               "gm90/analysis_hg38/BS-TAG-Sample-90_CpG_percentFix.methylKit",
               "gm91/analysis_hg38/BS-TAG-Sample-91_CpG_percentFix.methylKit",
               "gm92/analysis_hg38/BS-TAG-Sample-92_CpG_percentFix.methylKit",
               "gm93/analysis_hg38/BS-TAG-Sample-93_CpG_percentFix.methylKit")
output_dir ="Analysis_hg38/"

#' read in files, CpGs are unstranded and needs minimum of 10 coverage 
#+ read-files, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
familyObj=methRead(file.list, sample.id=list("gm77", "gm78", "gm79", "gm80", "gm81", 
"gm82", "gm83", "gm84", "gm85", "gm86", "gm87", "gm88", "gm89", "gm90", "gm91",
 "gm92", "gm93"), assembly="hg38_decoy_EBV", treatment=rep(0,17), context="CpG", mincov = 3)
 
#'  create methylation histograms
#+ Density-plot, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
for (i in 1:length(file.list)){
getMethylationStats(familyObj[[i]],plot=TRUE,both.strands=FALSE)
}
 
#' create coverage plots
#+ Coverage-plots, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
for (i in 1:length(file.list)){
getCoverageStats(familyObj[[i]],plot=TRUE,both.strands=FALSE)
}
 
#' destrand and require all samples to have the CpG covered
#+ unite, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
meth<-unite(object = familyObj, destrand=TRUE)
 
#' generate cluster dendogram
#+ cluster-plot, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
 
#' save cluster information for plotting with ggplot2 later
#+ cluster-save, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
saveRDS(hc,file="Analysis_hg38/dendogram_test2.rda")
 
#' plot of principal component contribution to variance 
#+ PC-plot, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
PCASamples(meth, screeplot=TRUE)

#' PCA plot
#+ PCA-plot, fig.width = 15, fig.height = 15, tidy.opts=list(blank=FALSE,width.cutoff=60)
PCASamples(meth)
 
#' read in and make lane and generation factors
#+ load-meta, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
metadata <- read.table("metadata-BS.txt", header=T)
metadata$Lane1 <- as.factor(metadata$Lane1)
metadata$generation <- as.factor(metadata$generation)
sampleAnnotation<- metadata[,c(3,4,5,7,8)]

#' perform association test between principal components and the metadata
#+ Association-analysis, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
as=assocComp(mBase=meth,sampleAnnotation)
 
#' Create a heatmap of the associations
#+ PCA-heatmap, include = TRUE, fig.width = 7, fig.height = 5, tidy.opts=list(blank=FALSE,width.cutoff=60)
hmcol3<-colorRampPalette(brewer.pal(9,"Blues"))(5)
hmcol4<-hmcol3
hmcol4[1]<-"white"
logpvals.raw <- -log10(as$association[,1:10])
heatmap.2(logpvals.raw,Rowv=F,Colv=rownames(logpvals.raw), dendrogram='none', trace='none',
          margins=c(8,8),colsep=c(1:11),rowsep=c(1:16),
          sepwidth=c(0.025,0.025), sepcolor="black", col=hmcol4,
		  breaks=c(0,1.30103,2,3,4,5),
          key.xlab = "-log10 pvalue", main="PC Assoc. heatmap")

#' Create some correlation plots
#+ Correlation-plots-all, include = TRUE, fig.width = 7, fig.height = 6, tidy.opts=list(blank=FALSE,width.cutoff=60)
mat <- percMethylation(meth)
cormat <- generate_heatmap(mat, "CpGs_test")
saveRDS(cormat, file=paste(output_dir,"cormat_all.rds", sep=""))

#' Filter the methylation file by canonical chromosomes as well as different coverage minimums
#+ filter-meth, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)  
idx_can <- grep("chr", meth$chr)
meth_can <- meth[idx_can,]
idx_cov_col <- grep(pattern = "coverage", x= names(meth))
idx_cov_col_Kids <- idx_cov_col[-c(1,2,13:16)]
meth_can_data <- getData(meth_can)

# Get number of CpGs passing various coverage thresholds for all 11 children
idx_min5 <- which(rowSums(meth_can_data[,idx_cov_col_Kids]>=5) == 11)
length(idx_min5)
idx_min10 <- which(rowSums(meth_can_data[,idx_cov_col_Kids]>=10) == 11)
length(idx_min10)
idx_min15 <- which(rowSums(meth_can_data[,idx_cov_col_Kids]>=15) == 11)
length(idx_min15)
idx_min20 <- which(rowSums(meth_can_data[,idx_cov_col_Kids]>=20) == 11)
length(idx_min20)

# save the meth databases at 10x, 15x, and 20x thresholds
meth_can_10 <- meth_can[idx_min10,]
saveRDS(meth_can_10, file=paste(output_dir,"meth_can_10.rds", sep=""))
meth_can_15 <- meth_can[idx_min15,]
saveRDS(meth_can_15, file=paste(output_dir,"meth_can_15.rds", sep=""))
meth_can_20 <- meth_can[idx_min20,]
saveRDS(meth_can_20, file=paste(output_dir,"meth_can_20.rds", sep=""))
saveRDS(meth,file=paste(output_dir,"meth_no_thres.rds", sep=""))

#' Create some correlation plots for different thresholds
#+ Correlation-plots-thresholds, include = TRUE, fig.width = 7, fig.height = 6, tidy.opts=list(blank=FALSE,width.cutoff=60)
mat_can_10<- percMethylation(meth_can_10)
cormat_can_10 <- generate_heatmap(mat_can_10, "threshold_10")
saveRDS(cormat_can_10, file=paste(output_dir,"cormat_can_10.rds", sep=""))
clusterSamples(meth_can_10, dist="correlation", method="ward", plot=TRUE)

mat_can_15<- percMethylation(meth_can_15)
cormat_can_15 <- generate_heatmap(mat_can_15, "threshold_15")
saveRDS(cormat_can_15, file=paste(output_dir,"cormat_can_15.rds", sep=""))
clusterSamples(meth_can_15, dist="correlation", method="ward", plot=TRUE)
 
mat_can_20<- percMethylation(meth_can_20)
cormat_can_20 <- generate_heatmap(mat_can_20, "threshold_20")
saveRDS(cormat_can_20, file=paste(output_dir,"cormat_can_20.rds", sep=""))
clusterSamples(meth_can_20, dist="correlation", method="ward", plot=TRUE)

#' Filter the previously filtered methylation file by sex chromosomes
#+ filter-sex, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)  
idx_X <- grep("chrX", meth_can$chr)
idx_Y <- grep("chrY", meth_can$chr)
idx_EBV <- grep("chrEBV", meth_can$chr)
meth_can_fil <- meth_can[-c(idx_X,idx_Y,idx_EBV),]
meth_can_fil_data <- getData(meth_can_fil)

# Get number of CpGs passing various coverage thresholds for all 11 children
idx_min5_fil <- which(rowSums(meth_can_fil_data[,idx_cov_col_Kids]>=5) == 11)
length(idx_min5_fil)
idx_min10_fil <- which(rowSums(meth_can_fil_data[,idx_cov_col_Kids]>=10) == 11)
length(idx_min10_fil)
idx_min15_fil <- which(rowSums(meth_can_fil_data[,idx_cov_col_Kids]>=15) == 11)
length(idx_min15_fil)
idx_min20_fil <- which(rowSums(meth_can_fil_data[,idx_cov_col_Kids]>=20) == 11)
length(idx_min20_fil)

# save the meth databases at 10x, 15x, and 20x thresholds
meth_can_fil_10 <- meth_can_fil[idx_min10_fil,]
saveRDS(meth_can_fil_10, file=paste(output_dir,"meth_can_fil_10.rds", sep=""))
meth_can_fil_15 <- meth_can_fil[idx_min15_fil,]
saveRDS(meth_can_fil_15, file=paste(output_dir,"meth_can_fil_15.rds", sep=""))
meth_can_fil_20 <- meth_can_fil[idx_min20_fil,]
saveRDS(meth_can_fil_20, file=paste(output_dir,"meth_can_fil_20.rds", sep=""))
saveRDS(meth_can_fil,file=paste(output_dir,"meth_can_fil.rds", sep=""))

#' Create some correlation plots for different thresholds
#+ Correlation-plots-fil-thresholds, include = TRUE, fig.width = 7, fig.height = 6, tidy.opts=list(blank=FALSE,width.cutoff=60)
mat_can_fil_10<- percMethylation(meth_can_fil_10)
cormat_can_10_fil <- generate_heatmap(mat_can_fil_10, "fil_threshold_10")
saveRDS(cormat_can_10_fil, file=paste(output_dir,"cormat_can_fil_10.rds", sep=""))

mat_can_15_fil <- percMethylation(meth_can_fil_15)
cormat_can_15_fil <- generate_heatmap(mat_can_15_fil, "fil_threshold_15")
saveRDS(cormat_can_15_fil, file=paste(output_dir,"cormat_can_fil_15.rds", sep=""))
 
mat_can_20_fil <- percMethylation(meth_can_fil_20)
cormat_can_20_fil <- generate_heatmap(mat_can_20_fil, "fil_threshold_20")
saveRDS(cormat_can_20_fil, file=paste(output_dir,"cormat_can_fil_20.rds", sep=""))



