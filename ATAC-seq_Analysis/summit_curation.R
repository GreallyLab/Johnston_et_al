# read in the merged summit files 
options(stringsAsFactors = F)
library(data.table)
summits <- fread("Summit_idrPeak_chr_250_noOver_sort_merge.txt", sep="\t")
summits <- as.data.frame(summits)
dim(summits) # 143,660
head(summits)
idx_summmit_NOoverlap <- grep(",", summits$V4, invert = T)
idx_summmit_overlap <- grep(",", summits$V4)

summit_NOoverlap <- summits[idx_summmit_NOoverlap,]
dim(summit_NOoverlap) # 13525
summit_overlap <- summits[idx_summmit_overlap,]
dim(summit_overlap) #130,135
str(summit_overlap)
i<-1
summary(summit_NOoverlap$stop-summit_NOoverlap$start)

sum((summit_overlap$V3-summit_overlap$V2)>1000 )
# 130135*2
# 260270+13525  276,590
260270+13525+2192*2
summary(summits$V3-summits$V2)
sum((summit_overlap$V3-summit_overlap$V2)>1000) # 2,192
which.max((summit_overlap$V3-summit_overlap$V2))
name_high_score <- rep(0, nrow(summit_overlap))
hist((summit_overlap$V3-summit_overlap$V2))
i<-1
i<-109
i<-59275
head(summit_overlap)
summit_overlap[i,]
new_summit_overlap<-NULL
1<-1
for (i in 1:nrow(summit_overlap)) {
#  vec_names <- strsplit(strsplit(summit_overlap$V4[[i]], ",")[[1]], "/")[[1]][3]
  vec_scores <- strsplit(summit_overlap$V5[[i]], ",")[[1]]
  if ((summit_overlap$V3[i]-summit_overlap$V2[i])<1000){
    new_start<- floor(mean(as.numeric(strsplit(summit_overlap$V6[[i]], ",")[[1]])))
    new_score <- mean(as.numeric(vec_scores))
    new_summit <- data.frame(chr=summit_overlap$V1[i], start=new_start-250, 
                             stop=new_start+250, name=summit_overlap$V4[i],
                             score=new_score, strand=".")
  }
  else {
    n <- floor((summit_overlap$V3[i]-summit_overlap$V2[i])/500)
    spacing <- ceiling((summit_overlap$V3[i]-summit_overlap$V2[i])/n)
    buffer <- floor(((summit_overlap$V3[i]-summit_overlap$V2[i])-(spacing*(n-1))-500)/2)
    new_start <- NULL
    for (j in 0:(n-1)) {
      new_start<- c(new_start, summit_overlap$V2[i]+buffer+(spacing*j))
    }
    new_score <- mean(as.numeric(vec_scores))
    new_summit <- data.frame(chr=summit_overlap$V1[i], start=new_start, 
                             stop=new_start+500, name=summit_overlap$V4[i],
                             score=new_score, strand=".")
  }
  new_summit_overlap <- rbind(new_summit_overlap, new_summit)
}
head(new_summit_overlap)
dim(new_summit_overlap) # 260,317
260317+13525 # 273842 !! need to be less than 276,590
head(new_summit_overlap)
table(new_summit_overlap$stop-new_summit_overlap$start)

colnames(summit_NOoverlap) <- colnames(new_summit_overlap)
cat_summit_NOoverlap <- rbind(summit_NOoverlap, new_summit_overlap)
write.table(cat_summit_NOoverlap, "Summit_idrPeak_chr_250_noOver_sort_merge_noOver.txt",
            append = F, quote = F, col.names = F, row.names = F, sep = "\t")

# fixing the .5s 
summit_bed <- fread("Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed")
summit_bed <- as.data.frame(summit_bed)
summit_bed$V2 <- floor(summit_bed$V2)
summit_bed$V3 <- floor(summit_bed$V3)

# fixing the 501 length summits 
idx_501 <- which((summit_bed$V3-summit_bed$V2)>500)
length(idx_501) #17958
summit_bed$V3[idx_501] <- summit_bed$V3[idx_501]-1
table((summit_bed$V3-summit_bed$V2))
write.table(summit_bed, file = "Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed", append = F,
            quote = F, sep = "\t", row.names = F, col.names = F)


# remove the non-canon chr and the sex chromosomes from the summit peaks!
summitPeak_bed <- fread("replicate_summit_peak_nonspecific_sort.bed")
summitPeak_bed <- as.data.frame(summitPeak_bed)
head(summitPeak_bed)
dim(summitPeak_bed) # 108115
# GL, KI, KN, JTFH, EBV, chrX, chrY
nonCanon <- c("GL", "KI", "KN", "JTFH", "EBV", "chrX", "chrY")
idx_rm_nonCanon_summitPeak<-NULL
for (i in 1:length(nonCanon)) {
  idx_rm_nonCanon_temp <- grep(pattern = nonCanon[i], 
                               x = summit_bed$V1) 
  idx_rm_nonCanon_summitPeak <- c(idx_rm_nonCanon_summitPeak,idx_rm_nonCanon_temp)
}
length(idx_rm_nonCanon_summitPeak) # 4424
head(idx_rm_nonCanon_summitPeak)
summitPeak_canon_bed <- summitPeak_bed[-idx_rm_nonCanon_summitPeak, ]
table(summitPeak_canon_bed$V1)
write.table(summitPeak_canon_bed, file = "replicate_summit_peak_nonspecific_sort_Canon.bed", append = F,
            quote = F, sep = "\t", row.names = F, col.names = F)


## still have overlaps?!
summit_filtered <- read.table("Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_merge.bed")
table(summit_filtered[,3]-summit_filtered[,2])
idx_dup <- which((summit_filtered[,3]-summit_filtered[,2])>500)
summit_filtered[idx_dup,4]

summit_merge <- read.table("Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack.bed")
head(summit_merge)
dim(summit_merge)

summit_dup_names <-NULL
for ( i in 1:length(idx_dup)){
  summit_dup_names <- c(summit_dup_names,strsplit(summit_filtered[idx_dup[i],4], ",")[[1]])
}
summit_merge_dup <- summit_merge[(summit_merge$V4 %in% summit_dup_names),]
head(summit_merge_dup)
summit_merge_dup$V2[1] <- summit_merge_dup$V2[1]-1
summit_merge_dup$V3[1] <- summit_merge_dup$V3[1]-1

summit_merge_dup$V2[3] <- summit_merge_dup$V2[3]-1
summit_merge_dup$V3[3] <- summit_merge_dup$V3[3]-1

summit_merge_dup$V2[5] <- summit_merge_dup$V2[5]-1
summit_merge_dup$V3[5] <- summit_merge_dup$V3[5]-1

summit_merge_dup$V2[7] <- summit_merge_dup$V2[7]+1
summit_merge_dup$V3[7] <- summit_merge_dup$V3[7]+1

summit_merge_dup$V2[8] <- summit_merge_dup$V2[8]-1
summit_merge_dup$V3[8] <- summit_merge_dup$V3[8]-1

summit_merge_dup$V2[10] <- summit_merge_dup$V2[10]-1
summit_merge_dup$V3[10] <- summit_merge_dup$V3[10]-1

summit_merge_dup$V2[12] <- summit_merge_dup$V2[12]-1
summit_merge_dup$V3[12] <- summit_merge_dup$V3[12]-1


summit_merge_nodup <- summit_merge[!(summit_merge$V4 %in% summit_dup_names),]
dim(summit_merge_nodup) #139575

summit_merge_fix <- rbind(summit_merge_dup,summit_merge_nodup)
dim(summit_merge_fix)
write.table(summit_merge_fix, "Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix.bed", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)

