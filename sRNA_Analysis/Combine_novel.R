novel_miRNA_bed <-read.table("../Merged_novel_fillist.bed", stringsAsFactors = F)
head(novel_miRNA_bed)
novel_miRNA_bed_reArrange <- novel_miRNA_bed[,c(1:3,6,5,4)]
head(novel_miRNA_bed_reArrange)
dim(novel_miRNA_bed_reArrange)
novel_miRNA_bed_reArrange_10 <- novel_miRNA_bed_reArrange[novel_miRNA_bed_reArrange$V5 > 10,]
head(novel_miRNA_bed_reArrange_10)
dim(novel_miRNA_bed_reArrange_10) # 78
write.table(novel_miRNA_bed_reArrange_10, "Individual_novel_miRNA.bed", append = F,
            quote = F, row.names = F, col.names = F, sep = "\t")

# combining the Individual and Pooled miRNAs 
novel_miRNA_bed_reArrange_10$Method <- "Individual"

overlaps <- read.table("Pooled_individual_overlap.txt")
head(overlaps)
idx_overlap_ind <- which(novel_miRNA_bed_reArrange_10$V6 %in% overlaps$V4)
novel_miRNA_bed_reArrange_10$Method[idx_overlap_ind] <- "Both"

# read in the pooled novel miRNA
Pooled <- read.table("Pooled_novel_miRNA.bed")
head(Pooled)
Pooled$Method <- "Pooled"
idx_overlap_pool <- which(Pooled$V4 %in% overlaps$V10)
Pooled_no_over <- Pooled[-idx_overlap_pool,]
head(Pooled_no_over)
dim(Pooled_no_over) #  143   7

colnames(novel_miRNA_bed_reArrange_10) <- colnames(Pooled_no_over)
All_novel_miRNA <- rbind(novel_miRNA_bed_reArrange_10, Pooled_no_over)
head(All_novel_miRNA)

# how many of the combined list share a name? 2 have the same name.. and I need to fix
table(table(All_novel_miRNA$V4))
# 1   2 
# 217   2 

# which two
which(duplicated(All_novel_miRNA$V4)) # 180 204 (both in the pooled, so make the name _pooled)

All_novel_miRNA[180,4] <- paste(All_novel_miRNA[180,4], "_pooled", sep = "")
All_novel_miRNA[204,4] <- paste(All_novel_miRNA[204,4], "_pooled", sep = "")

table(table(All_novel_miRNA$V4)) # all 221 now have different names

write.table(x = All_novel_miRNA[,-7], file = "Combined_novel_miRNA.bed", append = F, 
            quote = F, sep = "\t", row.names = F, col.names = F)

table(All_novel_miRNA$V1) # chrJTFH01000007.1, chrJTFH01000009.1 need to have the chr removed.





