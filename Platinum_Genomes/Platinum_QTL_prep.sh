# I want to make the prep files for calling QTLs with the new platinum data

# First I'll prepare the haplotype file
# removing the x chr from the haplotype file
# btw header for haplotype file is: 
#CHROM  START   STOP    NA12877 NA12878 NA12879 NA12880 NA12881 NA12882 NA12883 NA12884 NA12885 NA12886 NA12887 NA12888 NA12893
# also need to add the haplo number to each
awk '{print $0"\tHaplo_"NR-1}' haplotype_transmission_grch38.txt >  haplotype_transmission_grch38_named.txt
awk 'NR>1 && $1 != "chrX"' haplotype_transmission_grch38_named.txt | wc -l # 708
awk 'NR>1 && $1 != "chrX"' haplotype_transmission_grch38_named.txt > haplotype_transmission_grch38_named_noXYHeader.txt
# need to get a file that I can intersect with ATACseq peaks
awk '{print $1"\t"$2"\t"$3"\t"$17}' haplotype_transmission_grch38_named_noXYHeader.txt > haplotype_grch38_named_4intersect.txt

# bringing in ATAC peak file
cp ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_numPeak_noXY.bed .
awk '{print $1}' Rep_nonspecific_numPeak_noXY.bed | uniq -c 

awk '$1 ~ "chr1" || $1 ~ "chr2" || $1 == "chr3" || $1 == "chr4" || $1 == "chr5" || $1 == "chr6" || $1 == "chr7" || $1 == "chr8" || $1 == "chr9" {print $1}' Rep_nonspecific_numPeak_noXY.bed | wc -l # 101,741


module load bedtools2
bedtools intersect -a Rep_nonspecific_numPeak_noXY.bed -b haplotype_grch38_named_4intersect.txt -wo > ATAC_Plat_haplo_interesect.txt

wc -l ATAC_Plat_haplo_interesect.txt
# 101,537 ATAC_Plat_haplo_interesect.txt # ALMOST all peaks!

# How many previously intersecting with haplotype 
wc -l ../Data_Prev_Paper/Total_Peak_QTL/ATAC_haplo_intersect.txt 
# 93,089 ../Data_Prev_Paper/Total_Peak_QTL/ATAC_haplo_intersect.txt

bedtools intersect -a ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_numPeak_noXY.bed -b haplotype_grch38_named_4intersect.txt -wa | uniq -d
# chr12	21527129	21528377	Peak_19467	1000	.
# chr12	94872072	94874678	Peak_22522	1000	.
# chr19	1258829	1263587	Peak_40661	1000	.
# chr3	20047538	20050524	Peak_57154	1000	.


############### E-QTL prep ################
module load bedtools2
bedtools intersect -a gencode.v24.primary_assembly.annotation.EBV.GENES.bed -b haplotype_grch38_named_4intersect.txt -wa | uniq -d | wc -l 
# 357 deuplicated genes, so 357 that intersect two different haplotypes. 
# The genes on unplaced loci don't have chr in front of them and aren't included.. but ok because only looking at chr1-22

bedtools intersect -a gencode.v24.primary_assembly.annotation.EBV.GENES.bed -b haplotype_grch38_named_4intersect.txt -wo > RNA_haplo_intersect.txt
wc -l RNA_haplo_intersect.txt
# 57566 RNA_haplo_intersect.txt # so most of them intersected

/home/users/ajohnsto/17_member/RNA-seq/Kallisto/PCA_analysis/results_norm_filt_genes_noSex_rlog

############### splice-QTL prep ################
# For sQTL, we selected the most significant p-value among all transcripts for each gene. (Li et al.)

# need to generate a bed file of transcripts from the gtf file
awk 'BEGIN {OFS = "\t"}; $3=="transcript" {print $1"\t"$4"\t"$5"\t"$12"\t"$6"\t"$7}' gencode.v24.primary_assembly.annotation.EBV.gtf | wc -l
# 199476, which is correct

awk 'BEGIN {OFS = "\t"}; $3=="transcript" {print $1"\t"$4"\t"$5"\t"$12"\t"$6"\t"$7}' gencode.v24.primary_assembly.annotation.EBV.gtf | tr -d '";' > gencode.v24.primary_assembly.annotation.EBV.TRANS.bed
wc -l gencode.v24.primary_assembly.annotation.EBV.TRANS.bed
# 199,476 

module load bedtools2
bedtools intersect -a gencode.v24.primary_assembly.annotation.EBV.TRANS.bed -b haplotype_grch38_named_4intersect.txt -wa | uniq -d | wc -l 
# 1360 duplicated transcripts, so 1360 that intersect two different haplotypes. 
# The transcripts on unplaced loci don't have chr in front of them and aren't included.. but ok because only looking at chr1-22

bedtools intersect -a gencode.v24.primary_assembly.annotation.EBV.TRANS.bed -b haplotype_grch38_named_4intersect.txt -wo > Trans_haplo_intersect.txt
wc -l Trans_haplo_intersect.txt
# 193060 - most intersect 

##### sRNA (miRNA) QTL analysis
wc -l miRNA_count_final.bed
#707
module load bedtools2
bedtools intersect -a miRNA_count_final.bed -b haplotype_grch38_named_4intersect.txt -wa | uniq -d | wc -l # 0, which makes sense 

bedtools intersect -a miRNA_count_final.bed -b haplotype_grch38_named_4intersect.txt -wo > miRNA_haplo_intersect.txt
wc -l miRNA_haplo_intersect.txt
# 707 - all intersect 






