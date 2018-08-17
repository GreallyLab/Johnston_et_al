setwd("/Volumes/users-1/Andrew Johnston/17_member/ATACseq/Total_peak_analysis/")
peak_sum_peak_sum_idr_05 <- read.table("replicates_summit_peak_sum_idr_05")
peak_sum_idr_05 <- read.table("Reps_peaks_peak_sum_idr_05_chr")
head(peak_sum_peak_sum_idr_05)
dim(peak_sum_peak_sum_idr_05) # 108130 20

Diff_sumPeak <- log10(peak_sum_idr_05$V15/peak_sum_idr_05$V19)
Add_sumPeak <- peak_sum_idr_05$V15 + peak_sum_idr_05$V19
Avg_sumPeak <- log(Add_sumPeak/2)
length(Diff_sumPeak)
length(Add_sumPeak)
length(Avg_sumPeak)

library(ggplot2)

plot(Avg_sumPeak, Diff_sumPeak,ylab="log10(Rep1/Rep2)", xlab="ln((Rep1+Rep2)/2", ylim = c(-1,1))
abline(h=-.5, col="red")
abline(h=.5, col="red")

title("MAplot of Rep1 and Rep2")

Rep2_specific_sumPeak <- peak_sum_idr_05[Diff_sumPeak<(-.5),]
dim(Rep2_specific_sumPeak) #3
head(Rep2_specific_sumPeak)

write.table(Rep2_specific_sumPeak, file = "replicate_summit_idr05_Rep2specific_peaks", append = F,
            quote = F, sep = "\t", row.names = F, col.names = F)

Rep2_sumPeak_specific_peaks.bed <- Rep2_specific_sumPeak[,c(1:6)]
write.table(Rep2_sumPeak_specific_peaks.bed, file = "replicate_summit_idr05_Rep2specific_peaks.bed", append = F,
            quote = F, sep = "\t", row.names = F, col.names = F)

hist(Rep2_sumPeak_specific_peaks.bed$V5) # all are 1000 score

# Creating w/o peak list w/o rep2 specific peaks
Not_specific_sumPeak <- peak_sum_idr_05[Diff_sumPeak>(-.425),]
dim(Not_specific_sumPeak)
# [1] 108115     20

head(Not_specific_sumPeak)
Rep_nonspecific_peaks_sumPeak.bed <- Not_specific_sumPeak[,c(1:6)]
write.table(Rep_nonspecific_peaks_sumPeak.bed, file = "replicate_summit_peak_nonspecific.bed", append = F,
            quote = F, sep = "\t", row.names = F, col.names = F)
