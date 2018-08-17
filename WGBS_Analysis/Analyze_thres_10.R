#' This is a script to use the CpGs which pass a coverage >10x for all children in the family and do some preliminary QC.

#+ load-libs, include = TRUE
suppressPackageStartupMessages({
	library(methylKit)
	library(RColorBrewer)
	library(gplots)
	library(knitr)
})

#+ setup, include=FALSE
opts_chunk$set(fig.width = 7, fig.height = 7, fig.path='Figs_Thres_10/')

#+ setup-functions, include=FALSE
source("../Correlation_plots2.R")
options(scipen = 999)

#' Loading in the data
#+ load-in-data, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
meth_10_fil <- readRDS("meth_can_fil_10.rds")


#' still need to remove the y chromosome (error in script)
#+ remove-chrY, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
# just a check of chrX
idx_X <- grep("chrx", meth_10_fil$chr)
length(idx_X)
# removing the Y chr CpGs
idx_Y <- grep("chrY", meth_10_fil$chr)
length(idx_Y)
meth_10_fil <- meth_10_fil[-c(idx_Y),]

#' Then I removed the CpGs which were above the 99.9th percentile because they ould be potential PCR duplicates 
#+ remove-999, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
idx_cov_col <- grep(pattern = "coverage", x= names(meth_10_fil))
meth_10_fil_data <- getData(meth_10_fil)
perc_999 <- quantile(rowMeans(meth_10_fil_data[,idx_cov_col]), .999) #     99% = 117.9412
idx_999 <- which(rowMeans(meth_10_fil_data[,idx_cov_col])>perc_999)
length(idx_999)
meth_10_fil_no99 <- meth_10_fil[-idx_999,] 
nrow(meth_10_fil_no99)
saveRDS(meth_10_fil_no99, "meth_10_fil_no99.rds")

#' Getting the new correlation plots 
#+ Plot-cor-999, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60), fig.width = 7, fig.height = 6
mat_10_fil_no99<- percMethylation(meth_10_fil_no99)
nrow(mat_10_fil_no99)
cormat_can_10_fil_no99 <- generate_heatmap(mat_10_fil_no99, "fil_threshold_20")

#' Plotting the dendogram of the filtered CpGs
#+ Plot-dendogram-999, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
clusterSamples(meth_10_fil_no99, dist="correlation", method="ward", plot=TRUE)


#' Make a bed file for the filtered CpGs
#+ Make-bed, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
fil_10_99 <- data.frame(chr=meth_10_fil_no99$chr, start=meth_10_fil_no99$start, 
                        stop=meth_10_fil_no99$start+1, 
                        name=paste("CpG_",1:nrow(mat_10_fil_no99), sep=""),
                        score=".", strand=".") 
dim(fil_10_99)
write.table(fil_10_99, file="fil_10_99.bed",append = F, quote = F, 
            row.names = F,col.names = F, sep = "\t")
saveRDS(fil_10_99, "fil_10_99.rds")
			
#' Plotting the methylation density
#+ Plot-meth-denisty, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60), fig.width = 8, fig.height = 5
meth_t <- data.frame(IDs=colnames(mat_10_fil_no99), t(mat_10_fil_no99))
meth_t_melt <- melt(meth_t, id = "IDs")
#pdf("MethPerc_density_fil_20_99.pdf", width = 8, height=5)
ggplot(data = meth_t_melt, aes(x = value, col = IDs)) + 
	geom_density(aes(size = IDs)) + 
	ggtitle(paste("Methylation distribution-not smoothed-", "20-fil-99" , sep="")) + 
	xlab("CpG Methylation") + 
	ylab("Density") + theme_bw() +
	scale_size_manual(values = c(rep(.5, 12), .5, 1,1, rep(.5, 2)))


#' Plotting the average chromosomal coverage for each individual
#+ Plot-chr-avg-coverage, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60), fig.width = 10, fig.height = 6
meth_10_fil_no99_data <- getData(meth_10_fil_no99)
meth_10_fil_no99_data_spl <- split(meth_10_fil_no99_data, meth_10_fil_no99_data$chr)
meth_10_fil_no99_coverage <- sapply(meth_10_fil_no99_data_spl, function(x) {
  cov <- colMeans(x[,idx_cov_col])
  names(cov) <- colnames(mat_10_fil_no99)
  cov
})
meth_10_fil_no99_coverage <- as.data.frame(meth_10_fil_no99_coverage)
meth_10_fil_no99_coverage$Samples <- rownames(meth_10_fil_no99_coverage)
summary(meth_10_fil_no99_coverage$chr21)
melt_coverage <- melt(meth_10_fil_no99_coverage, id = "Samples")
ggplot(data = melt_coverage, aes(Samples, value)) + geom_boxplot() +
  ggtitle("Average Chromosomal Coverage by Individual") + 
  ylab("Average CpG Coverage") + 
  xlab("Individual")
ggplot(data = melt_coverage, aes(variable, value)) + geom_boxplot() +
  geom_jitter(aes(col=Samples)) +
  ggtitle("Average Chromosomal Coverage by Chromosome") + 
  ylab("Average CpG Coverage") + 
  xlab("Chromosome") +
  scale_colour_discrete(name="Individual")

#' get the intermediate methylation sites 
#+ Make-IM-bed, include = TRUE, tidy.opts=list(blank=FALSE,width.cutoff=60)
meth_10_fil_no99_data <- getData(meth_10_fil_no99)
mat_10_fil_no99<- percMethylation(meth_10_fil_no99)
idx_meth_col_Kids <- c(3:12,7)
idx_IM_kids <-  which(rowSums(mat_10_fil_no99[,idx_meth_col_Kids] <= 25 &
                                mat_10_fil_no99[,idx_meth_col_Kids] <= 75) >= 3)
length(idx_IM_kids)
meth_10_fil_no99_IM <- meth_10_fil_no99[idx_IM_kids,]

meth_10_fil_no99_IM_bed <- data.frame(chr=meth_10_fil_no99_IM$chr, start=meth_10_fil_no99_IM$start, 
                        stop=meth_10_fil_no99_IM$start+1, 
                        name=paste("CpG_",1:nrow(meth_10_fil_no99_IM), sep=""),
                        score=".", strand=".") 
dim(meth_10_fil_no99_IM_bed)
write.table(meth_10_fil_no99_IM_bed, file="fil_10_99_IM.bed",append = F, quote = F, 
            row.names = F,col.names = F, sep = "\t")
saveRDS(meth_10_fil_no99_IM_bed, "fil_10_99_IM.rds")
	
	

