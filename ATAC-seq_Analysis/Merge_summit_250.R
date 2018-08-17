setwd("/Volumes/users-1/Andrew Johnston/17_member/ATACseq/Rep1/Fastq/comb_reads/Peaks_shift/")
# read in the merged summit files 
options(stringsAsFactors = F)
library(data.table)
summits <- fread("Rep1_summits_chr_250_merge.bed", sep="\t")
summits <- as.data.frame(summits)
dim(summits) # 233,314
head(summits)
idx_summmit_NOoverlap <- grep(",", summits$V4, invert = T)
idx_summmit_overlap <- grep(",", summits$V4)


summit_NOoverlap <- summits[idx_summmit_NOoverlap,]
dim(summit_NOoverlap) # 198756
summit_overlap <- summits[idx_summmit_overlap,]
dim(summit_overlap) #34,558
str(summit_overlap)
i<-1
name_high_score <- rep(0, nrow(summit_overlap))
for (i in 1:nrow(summit_overlap)) {
  vec_names <- strsplit(summit_overlap$V4[[i]], ",")[[1]]
  vec_scores <- strsplit(summit_overlap$V5[[i]], ",")[[1]]
  name_high_score[i] <- vec_names[which.max(vec_scores)]
}
sum(name_high_score==0)
head(name_high_score)
length(name_high_score) # 34558

# read in all the 250 peaks,
all_summit <- fread("Rep1_summits_chr_250.bed", sep = "\t")
all_summit <- as.data.frame(all_summit)
head(all_summit)
dim(all_summit) # 277,900
all_summit_over_highScore <- all_summit[(all_summit$V4 %in% name_high_score),]
dim(all_summit_over_highScore) # 34558

cat_summit_noOver <- rbind(summit_NOoverlap,all_summit_over_highScore)
dim(cat_summit_noOver) # 233,314

write.table(x = cat_summit_noOver, file = "Rep1_summits_chr_250_noOver.bed", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)
