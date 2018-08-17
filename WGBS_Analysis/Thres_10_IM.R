meth_10_fil_no99 <- readRDS("meth_10_fil_no99.rds")
data_meth_10_fil_no99 <- getData(meth_10_fil_no99)
head(data_meth_10_fil_no99)
mat_meth_10_fil_no99 <- percMethylation(meth_10_fil_no99)
head(mat_meth_10_fil_no99)
idx_meth_col_Kids <- c(3:12,7)
idx_IM_kids <-  which(rowSums(mat_meth_10_fil_no99[,idx_meth_col_Kids] >= 25 &
                                mat_meth_10_fil_no99[,idx_meth_col_Kids] <= 75) >= 3)
length(idx_IM_kids) #5,395,058
dim(fil_10_99_IM) # 5,395,058 


meth_10_fil_no99_IM <- meth_10_fil_no99[idx_IM_kids,]
saveRDS(meth_10_fil_no99_IM, "meth_10_fil_no99_IM.rds")
meth_10_fil_no99_IM
meth_10_fil_no99_IM_bed <- data.frame(chr=meth_10_fil_no99_IM$chr, start=meth_10_fil_no99_IM$start, 
                                      stop=meth_10_fil_no99_IM$start+1, 
                                      name=paste("CpG_",1:nrow(meth_10_fil_no99_IM), sep=""),
                                      score=".", strand=".") 

write.table(meth_10_fil_no99_IM_bed, file="fil_10_99_IM.bed",append = F, quote = F, 
            row.names = F,col.names = F, sep = "\t")
saveRDS(meth_10_fil_no99_IM_bed, "meth_10_fil_no99_IM_bed.rds")
head(meth_10_fil_no99_IM_bed)
dim(meth_10_fil_no99_IM_bed)

mat_meth_10_fil_no99_IM <- mat_meth_10_fil_no99[idx_IM_kids,]
dim(mat_meth_10_fil_no99_IM) # 5,395,058
head(mat_meth_10_fil_no99_IM)
head(meth_10_fil_no99_IM)
meth_20_fil_no99_IM_meth_perc <- cbind(meth_10_fil_no99_IM_bed,mat_meth_10_fil_no99_IM)
dim(meth_20_fil_no99_IM_meth_perc) # 5395058 23
chrQTL_peaks_IM <- read.table("../../QTLs_2017/Integration/chrQTL_peaks_IM.txt")
head(chrQTL_peaks_IM)
dim(chrQTL_peaks_IM)


meth_20_fil_no99_IM_meth_perc_peak <- 
  meth_20_fil_no99_IM_meth_perc[meth_20_fil_no99_IM_meth_perc$name %in% chrQTL_peaks_IM$V10,]

meth_20_fil_no99_IM_meth_perc_peak <-
  merge(chrQTL_peaks_IM, meth_20_fil_no99_IM_meth_perc, by.x="V10", by.y="name")
head(meth_20_fil_no99_IM_meth_perc_peak)
dim(meth_20_fil_no99_IM_meth_perc_peak)
saveRDS(meth_20_fil_no99_IM_meth_perc_peak,"meth_20_fil_no99_IM_meth_perc_peak.rds")
head(meth_20_fil_no99_IM)


#### Need to get the CpGs underneath the summits 

meth_10_fil_no99_IM <- readRDS( "meth_10_fil_no99_IM.rds")
meth_10_fil_no99_IM <- data.frame(chr=meth_10_fil_no99_IM$chr, start=meth_10_fil_no99_IM$start, 
                                  stop=meth_10_fil_no99_IM$start+1, 
                                  name=paste("CpG_",1:nrow(meth_10_fil_no99_IM), sep=""),
                                  score=".", strand=".") 
mat_meth_10_fil_no99_IM <- percMethylation([meth_10_fil_no99_IM])
chrQTL_summits_IM <- read.table("../../QTLs_2017/Integration/Summits_intersect_HetVars_IM.txt")
head(chrQTL_summits_IM)
dim(chrQTL_summits_IM) # 362  13

head(meth_20_fil_no99_IM_meth_perc)

meth_20_fil_no99_IM_meth_perc_summit <-
  merge(chrQTL_summits_IM, meth_20_fil_no99_IM_meth_perc, by.x="V10", by.y="name")
head(meth_20_fil_no99_IM_meth_perc_summit)
dim(meth_20_fil_no99_IM_meth_perc_summit)
saveRDS(meth_20_fil_no99_IM_meth_perc_summit,"meth_20_fil_no99_IM_meth_perc_summit.rds")


