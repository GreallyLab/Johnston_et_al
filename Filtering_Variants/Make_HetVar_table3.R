setwd("/Volumes/home-1/anjohnst/17_member/QTLs_2017/Integration")
chrQTL <- readRDS("../ATAC_summit_QTL_2/ATAC_haplo_chrQTL_edit_matpat.rds")
#chrQTL <- readRDS("../ATAC_summit_QTL_2/Child_summit_5.bed")
options(stringsAsFactors = F)


head(chrQTL)

# Creating the HetVar in chrQTL peak table
 # this is current 
HetVar_in_peaks<- read.table("HetVar_intersect_chrQTL_peaks.txt")
head(HetVar_in_peaks)
dim(HetVar_in_peaks) #923
HetVar_table <- HetVar_in_peaks[,c(1:6,10)]
head(HetVar_table)
# include the peak Intermediate methylation information
peak_IM <- read.table("chrQTL_peaks_IM.txt")
head(peak_IM)
dim(peak_IM)
str(peak_IM)

HetVar_table$Peak_IM <- 0
HetVar_table$Peak_IM[HetVar_table$V10 %in% peak_IM$V4] <- 1
table(HetVar_table$Peak_IM) # 762

#HetVar in Summits 
summit_HV <- read.table("HetVars_intersect_summit.txt")
head(summit_HV)
dim(summit_HV)
str(summit_HV)
head(HetVar_table)
class(HetVar_table$V4)
HetVar_table$V4 <- as.character(HetVar_table$V4)
HetVar_table$Summit <- 0
HetVar_table$Summit[HetVar_table$V4 %in% summit_HV$V4] <- 1
table(HetVar_table$Summit) # 321
idx_summit <- which(summit_HV$V4 %in% HetVar_table$V4)
length(idx_summit)
summit_HV_only <- summit_HV[-idx_summit,]
dim(summit_HV_only)
summit_HV_only_rdy <- data.frame(summit_HV_only[,c(1:6)],V10="Summit_only", Peak_IM=0, Summit=1)
head(HetVar_table)
HetVar_table <- rbind(HetVar_table, summit_HV_only_rdy)
dim(HetVar_table)

# HetVar Summit_name
summit_HV_in_peak <- summit_HV[summit_HV$V4 %in% HetVar_table$V4,]
head(summit_HV_in_peak)
dim(summit_HV_in_peak)
summit_HV_in_peak_order <- match(HetVar_table$V4[HetVar_table$V4 %in% summit_HV$V4],summit_HV_in_peak$V4)
tail(summit_HV_in_peak_order)
HetVar_table$Summit_name <- "None"
HetVar_table$Summit_name[HetVar_table$V4 %in% summit_HV$V4] <- as.character(summit_HV_in_peak$V10[summit_HV_in_peak_order])
tail(HetVar_table,10)
dim(HetVar_table)

# HetVar Summit IM
summit_HV_IM <- read.table("Summits_intersect_HetVars_IM.txt")
head(summit_HV_IM)
dim(summit_HV_IM)
HetVar_table$Summit_IM <- 0
HetVar_table$Summit_IM[HetVar_table$Summit_name %in% summit_HV_IM$V4] <- 1
table(HetVar_table$Summit_IM) #138
dim(HetVar_table)

# HetVar eQTL
# first need to get the HetVar haplo information
HetVar_bed <- HetVar_table[,c(1:6)]
head(HetVar_bed)
write.table(HetVar_bed, "HetVar_in_chrQTL_peaks_2.bed", append = F, quote = F,
            col.names = F, row.names = F, sep = "\t")
#bedtools intersect -a HetVar_in_chrQTL_peaks_2.bed -b ../haplotype_grch38_named_4intersect.txt -wo > HetVar_in_chrQTL_peaks_haplo2.txt
HetVar_haplo <- read.table("HetVar_in_chrQTL_peaks_haplo2.txt")
head(HetVar_haplo)
HetVar_haplo[HetVar_haplo=="Haplo_31"]
dim(HetVar_haplo)
HetVar_table$Haplo <- HetVar_haplo$V10
head(HetVar_table)
tail(HetVar_table)
dim(HetVar_table)

# Adding in associated gene-eQTLs and miRNA eQTLs 
Gene_eQTLs<- readRDS("../../RNA-seq/Quant_2018/Genes_haplo_protlinc_sig2.rds")
head(Gene_eQTLs)

Gene_eQTLs[Gene_eQTLs$Gene_sym=="GSTM1",]

HetVar_table[HetVar_table$Haplo=="Haplo_31",]

HetVar_table$Gene_eQTL <- 0
HetVar_table$Gene_eQTL[HetVar_table$Haplo %in% Gene_eQTLs$V10] <- 1
table(HetVar_table$Gene_eQTL) # 479

# miRNA 
miRNA_eQTL <- readRDS("../miRNA_QTLs/miRNA_haplo_eQTL_edit_matpat.rds")
head(miRNA_eQTL)
dim(miRNA_eQTL) # 47  25

HetVar_table$miRNA_eQTL <- 0
HetVar_table$miRNA_eQTL[HetVar_table$Haplo %in% miRNA_eQTL$V10] <- 1
table(HetVar_table$miRNA_eQTL) #119

# adding in the gene names for each variant
Gene_eQTLs_short <- Gene_eQTLs[,c(2,1)]
table(table(Gene_eQTLs_short$V10)) # 6 on 1

spl_gene_eQTL_short <- split(Gene_eQTLs_short, f = Gene_eQTLs_short$V10)
genes_in_haplo <- sapply(spl_gene_eQTL_short, FUN = function(x) {
  if(nrow(x)>1) {
    temp <- x[1,1]
    out <- data.frame(haplo=x[1,2], genes=paste(temp,x[2,1], sep=","))
  }
  else {
    out <- data.frame(haplo=x[1,2], genes=x[1,1])
  }
})
t_genes_in_haplo <- t(genes_in_haplo)
head(t_genes_in_haplo)
dim(t_genes_in_haplo) #194
t_genes_in_haplo <- as.data.frame(t_genes_in_haplo)
dim(t_genes_in_haplo) #194 
rownames(t_genes_in_haplo) <- NULL
head(HetVar_table)
HetVar_table_merge_gene <- merge(HetVar_table, t_genes_in_haplo, by.x ="Haplo",
                                 by.y="haplo", all.x = T)
head(HetVar_table_merge_gene)
HetVar_table_merge_gene$genes <- unlist(HetVar_table_merge_gene$genes)
class(HetVar_table_merge_gene$genes)
sum(!is.na(HetVar_table_merge_gene$genes)) #479 variants
HetVar_table_merge_gene$genes <- replace(x = HetVar_table_merge_gene$genes,
                                         list = is.na(HetVar_table_merge_gene$genes), 
                                         values = "None")
class(HetVar_table_merge_gene$genes)
table((HetVar_table_merge_gene$genes))
str(HetVar_table_merge_gene)

head(miRNA_eQTL)
table(table(miRNA_eQTL$V10)) # 39 + 4 
mi_eQTL_short <- miRNA_eQTL[,c(4,7)]
table(table(mi_eQTL_short$V10)) # 39 + 4 (so 4 on two)

spl_mi_eQTL_short <- split(mi_eQTL_short, f = mi_eQTL_short$V10)
miRNA_in_haplo <- sapply(spl_mi_eQTL_short, FUN = function(x) {
  if(nrow(x)>1) {
    temp <- x[1,1]
    out <- data.frame(haplo=x[1,2], miRNA=paste(temp,x[2,1], sep=","))
  }
  else {
    out <- data.frame(haplo=x[1,2], miRNA=x[1,1])
  }
})
t_miRNA_in_haplo <- t(miRNA_in_haplo)
head(t_miRNA_in_haplo)

HetVar_table_merge_gene_miRNA <- merge(HetVar_table_merge_gene, t_miRNA_in_haplo, by.x ="Haplo",
                                 by.y="haplo", all.x = T)
head(HetVar_table_merge_gene_miRNA)
HetVar_table_merge_gene_miRNA$miRNA <- unlist(HetVar_table_merge_gene_miRNA$miRNA)
sum(!is.na(HetVar_table_merge_gene_miRNA$miRNA)) #119 variants
HetVar_table_merge_gene_miRNA$miRNA <- replace(x = HetVar_table_merge_gene_miRNA$miRNA,
                                         list = is.na(HetVar_table_merge_gene_miRNA$miRNA), 
                                         values = "None")
class(HetVar_table_merge_gene_miRNA$miRNA)
table((HetVar_table_merge_gene_miRNA$miRNA))

# Finding out if the HetVar's check out according to matPat discovery of chrQTL
dim(HetVar_table_merge_gene_miRNA) # 931  16
head(HetVar_table_merge_gene_miRNA)
saveRDS(HetVar_table_merge_gene_miRNA, "HetVar_table_merge_gene_miRNA3.rds")
# HetVar_table_merge_gene_miRNA <- readRDS("HetVar_table_merge_gene_miRNA2.rds")
# loading in the vcf file for Heterozygous variants
library(data.table)
All_HetVar_vcf <- fread("../../Platinum_Genomes/files/GenotypeFiles/phg000812.v1.PlatinumGenomes.genotype-calls-vcf.c1/High_Confidence_Calls_HG38_hetero.vcf", sep="\t")
head(All_HetVar_vcf)
dim(All_HetVar_vcf) # 4216912      22

All_HetVar_vcf$Name <- paste("Variant_",seq(1,4216912,1),sep = "")
head(All_HetVar_vcf)

HetVar_in_peaks_merge_geno <- merge(HetVar_table_merge_gene_miRNA, All_HetVar_vcf, 
                                    by.x = "V4", by.y = "Name", all.x = T)
head(HetVar_in_peaks_merge_geno)
dim(HetVar_in_peaks_merge_geno)
HetVar_in_peaks_merge_geno_short <- HetVar_in_peaks_merge_geno[,c(1,20:21,26:38)]
head(HetVar_in_peaks_merge_geno_short) 
dim(HetVar_in_peaks_merge_geno_short) #  931  16

#### merge summit allele with peaks
df_sig_ATAC_names <- readRDS("../ATAC_summit_QTL_2/df_sig_ATAC_names.rds")
head(df_sig_ATAC_names)
Summit_chrQTL_intersect_peaks <- read.table("Summit_chrQTL_intersect_peaks.txt")
head(Summit_chrQTL_intersect_peaks)
Summit_chrQTL_intersect_peaks_short <- Summit_chrQTL_intersect_peaks[,c(4,10)]
sum(duplicated(Summit_chrQTL_intersect_peaks_short$V10)) #7
head(Summit_chrQTL_intersect_peaks_short)
Summit_chrQTL_intersect_peaks_short_merge_sig_names <- 
  merge(Summit_chrQTL_intersect_peaks_short, df_sig_ATAC_names, by.x = "V4",
        by.y = "Summit")
head(Summit_chrQTL_intersect_peaks_short_merge_sig_names)
dim(Summit_chrQTL_intersect_peaks_short_merge_sig_names) # 276  3

Summit_chrQTL_intersect_peaks_short_merge_sig_names
peak_ATAC_sig_name <- Summit_chrQTL_intersect_peaks_short_merge_sig_names[,c(2,3)]
head(peak_ATAC_sig_name)
dim(peak_ATAC_sig_name) # 276
sum(duplicated(peak_ATAC_sig_name$V10)) #7
idx_dup_sigName <- duplicated(peak_ATAC_sig_name$V10)
dup_sigName <- peak_ATAC_sig_name$V10[idx_dup_sigName]
peak_ATAC_sig_name_noDup <- peak_ATAC_sig_name[!idx_dup_sigName,]
dim(peak_ATAC_sig_name_noDup) # 269
i<-1
for (i in 1:length(dup_sigName)) {
  dups <- peak_ATAC_sig_name[peak_ATAC_sig_name$V10 == dup_sigName[i],]
  if(dups[1,2] == "Both" | dups[2,2] == "Both") {
    peak_ATAC_sig_name_noDup[peak_ATAC_sig_name_noDup$V10 == dup_sigName[i],2] <- "Both"
  }
  if(dups[1,2] == "Mat" & dups[2,2] == "Mat") {
    peak_ATAC_sig_name_noDup[peak_ATAC_sig_name_noDup$V10 == dup_sigName[i],2] <- "Mat"
  }
  if(dups[1,2] == "Mat" & dups[2,2] == "Pat") {
    peak_ATAC_sig_name_noDup[peak_ATAC_sig_name_noDup$V10 == dup_sigName[i],2] <- "Both"
  }
  if(dups[1,2] == "Pat" & dups[2,2] == "Pat") {
    peak_ATAC_sig_name_noDup[peak_ATAC_sig_name_noDup$V10 == dup_sigName[i],2] <- "Pat"
  }
  if(dups[1,2] == "Pat" & dups[2,2] == "Mat") {
    peak_ATAC_sig_name_noDup[peak_ATAC_sig_name_noDup$V10 == dup_sigName[i],2] <- "Both"
  }
}
peak_ATAC_sig_name_noDup[peak_ATAC_sig_name_noDup$V10 %in% dup_sigName,]

head(HetVar_table_merge_gene_miRNA)

HetVar_table_merge_gene_miRNA[HetVar_table_merge_gene_miRNA$V10=="Summit_only",]
Summit_chrQTL_intersect_peaks_short_merge_sig_names[Summit_chrQTL_intersect_peaks_short_merge_sig_names$V4 %in% summit_HV_only$V10,]

HetVar_table_merge_gene_miRNA_summitonly<- HetVar_table_merge_gene_miRNA[HetVar_table_merge_gene_miRNA$V10=="Summit_only",]

HetVar_table_merge_gene_miRNA_sigName <- merge(HetVar_table_merge_gene_miRNA, 
                                               peak_ATAC_sig_name_noDup, by = "V10", all.x = T)
head(HetVar_table_merge_gene_miRNA_sigName)                          
dim(HetVar_table_merge_gene_miRNA_sigName) 
order_allele <- match(HetVar_table_merge_gene_miRNA_sigName$Summit_name[is.na(HetVar_table_merge_gene_miRNA_sigName$Allele)],
                      Summit_chrQTL_intersect_peaks_short_merge_sig_names$V4)
HetVar_table_merge_gene_miRNA_sigName$Allele[is.na(HetVar_table_merge_gene_miRNA_sigName$Allele)] <- Summit_chrQTL_intersect_peaks_short_merge_sig_names$Allele[order_allele]
tail(HetVar_table_merge_gene_miRNA_sigName)

HetVar_table_merge_gene_miRNA_sigName_merge_Haplo <- 
  merge(HetVar_table_merge_gene_miRNA_sigName, HetVar_in_peaks_merge_geno_short, by="V4")

head(HetVar_table_merge_gene_miRNA_sigName_merge_Haplo)
dim(HetVar_table_merge_gene_miRNA_sigName_merge_Haplo) # 931  32
HetVar_table_final <- HetVar_table_merge_gene_miRNA_sigName_merge_Haplo[,c(4:6,1,7:8,2:3,9:32)]
head(HetVar_table_final)

HetVar_table_final$chrQTL_check <- 0
for(i in 1:nrow(HetVar_table_final)) {
  if (HetVar_table_final$Allele[i] == "Pat" & 
      (HetVar_table_final$NA12877[i] == "1|0" | HetVar_table_final$NA12877[i] == "0|1")) {
    HetVar_table_final$chrQTL_check[i] <- 1
  }
  if (HetVar_table_final$Allele[i] == "Mat" & 
      (HetVar_table_final$NA12878[i] == "1|0" | HetVar_table_final$NA12878[i] == "0|1")) {
    HetVar_table_final$chrQTL_check[i] <- 1
  }
  if (HetVar_table_final$Allele[i] == "Both" & 
      (HetVar_table_final$NA12877[i] == "1|0" | HetVar_table_final$NA12877[i] == "0|1") &
      (HetVar_table_final$NA12878[i] == "1|0" | HetVar_table_final$NA12878[i] == "0|1")) {
    HetVar_table_final$chrQTL_check[i] <- 1
  }
}
table(HetVar_table_final$chrQTL_check) # 313 618 
head(HetVar_table_final[HetVar_table_final$Allele=="Both",])
saveRDS(HetVar_table_final,"HetVar_table_final3.rds")
#HetVar_table_final <- readRDS("HetVar_table_final3.rds")
head(HetVar_table_final)

HetVar_table_final[HetVar_table_final$genes=="ENSG00000134184.12",]



## adding in the Li-et-al overlaps
li_ID_overlap <- read.table("../../RNA-seq/Quant_2018/sig_qvals_merge_protlinc_overLi2.txt", header=T)
li_ID_overlap <- read.table("../../RNA-seq/Quant_2018/sig_qvals_merge_protlinc_overLi_IDonly2.txt", header=T)

head(li_ID_overlap)
dim(li_ID_overlap)

Gene_list_li <- data.frame(ID=as.numeric(), gene=as.character())
i<-21
for(i in 1:nrow(HetVar_table_final)){
    temp_genes <- HetVar_table_final$genes[i]
    temp_genes_c <- do.call(rbind, strsplit(temp_genes,","))[1,]
    temp_genes_df <- data.frame(ID=i, gene=temp_genes_c)
    Gene_list_li<-rbind(Gene_list_li, temp_genes_df)
}
     
dim(Gene_list_li) # 1240
head(Gene_list_li)
tail(Gene_list_li)
sum(Gene_list_li$gene %in% li_ID_overlap$Genes) #51 vs 48
idx_li_Gover <- which(Gene_list_li$gene %in% li_ID_overlap$Genes) #51
table(table(Gene_list_li[idx_li_Gover,1])) # only ever one gene on a haplotype block
Gene_list_li[idx_li_Gover,] # only ever one gene on a haplotype block


idx_ID_over <- Gene_list_li[idx_li_Gover,1]
HetVar_table_final$Li_overlap <- 0
HetVar_table_final$Li_overlap[idx_ID_over] <- 1

# analyzing overlaps

sum(HetVar_table_final$chrQTL_check==1) # 618

sum(HetVar_table_final$Gene_eQTL==1 &
      HetVar_table_final$chrQTL_check==1) #311

sum(HetVar_table_final$miRNA_eQTL==1 &
      HetVar_table_final$chrQTL_check==1) #86

sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1) #311
sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      (HetVar_table_final$Peak_IM == 1 | HetVar_table_final$Summit_IM == 1)) #253

sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      (HetVar_table_final$Peak_IM == 1 & HetVar_table_final$Summit_IM == 1)) #261

sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Peak_IM == 1 & HetVar_table_final$Li_overlap==1) #28

idx_li_varOver <- which(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Peak_IM == 1 & HetVar_table_final$Li_overlap==1) #28
HetVar_table_final[idx_li_varOver,15]

geneQTL_peakIM_li_et_al <- HetVar_table_final[idx_li_varOver,c(1:4)]
write.table(geneQTL_peakIM_li_et_al, "geneQTL_peakIM_li_et_al2.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
test <- read.table("geneQTL_peakIM_li_et_al2.txt")
test
# only variants from 8 genes are implicated: 
#GBP3 (no summit;GTF2B), SMN2(NAIP), DFNA5/GSDME, [ENTPD1, ACCS (no summit), not found when only ENSIDs filtering], RPS26, TBC1D4, POLI. Each gene has at least one variant that resides within the implicated summit except GBP3 and ACCS. 
# the variants that we tested were within the summit regions of TBC1D4.
# There was no summit intermediate methylation but I don't believe any CpGs fall within the summit region. 

# how many rare variants are present?
options(stringsAsFactors = FALSE)
rare_vars <- read.table("rare_variants.txt")
str(rare_vars)
str(geneQTL_peakIM_li_et_al)
sum(rare_vars$V1 %in% geneQTL_peakIM_li_et_al$V4) #1 
rare_vars[rare_vars$V1 %in% geneQTL_peakIM_li_et_al$V4,]
HetVar_table_final[idx_li_varOver,][(HetVar_table_final[idx_li_varOver,4] %in% rare_vars$V1),]




 sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Peak_IM == 1 & HetVar_table_final$Summit_IM == 0) #211
sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Summit_IM == 1) #50

idx_geneQTL_PeakIM <- which(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Peak_IM == 1 & HetVar_table_final$Summit_IM == 0) #202
idx_geneQTL_sumIM <- which(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Summit_IM == 1) #51

geneQTL_PeakIM_vars <- HetVar_table_final[idx_geneQTL_PeakIM,c(1:4)]
geneQTL_sumIM_vars <- HetVar_table_final[idx_geneQTL_sumIM,c(1:4)]
write.table(geneQTL_PeakIM_vars, "eQTL_peakIM_check_vars2.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
write.table(geneQTL_sumIM_vars, "eQTL_sumIM_check_vars2.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
check_sum_1 <- read.table("eQTL_sumIM_check_vars2.txt")
head(check_sum_1)
sum(check_sum_1$V4 %in% geneQTL_sumIM_vars$V4) #49

check_sum_1[!(check_sum_1$V4 %in% geneQTL_sumIM_vars$V4),]
geneQTL_sumIM_vars[!(geneQTL_sumIM_vars$V4 %in% check_sum_1$V4),]

check_peak_1 <- read.table("eQTL_peakIM_check_vars2.txt")
sum(check_peak_1$V4 %in% geneQTL_PeakIM_vars$V4) #202
check_peak_1[!(check_peak_1$V4 %in% geneQTL_PeakIM_vars$V4),]
geneQTL_PeakIM_vars[!(geneQTL_PeakIM_vars$V4 %in% check_peak_1$V4),]


#############################################



SNPtable <- fread("SNP150.txt")




sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Peak_IM == 1 & HetVar_table_final$Summit_IM == 0 &
      HetVar_table_final$miRNA_eQTL==0) #144
sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$chrQTL_check==1 &
      HetVar_table_final$Summit_IM == 1 & HetVar_table_final$miRNA_eQTL==0) #33


sum((HetVar_table_final$miRNA_eQTL==1 |HetVar_table_final$miRNA_eQTL==1) &
      HetVar_table_final$chrQTL_check==1 & HetVar_table_final$Peak_IM == 1) #79



# Loci that I previously checked:
sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$Summit_IM == 1 & 
  HetVar_table_final$chrQTL_check==1 & HetVar_table_final$miRNA_eQTL == 0)

sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$Summit_IM == 1 & 
      HetVar_table_final$chrQTL_check==1) #10
idx_eQTL_sumIM_check <- which(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$Summit_IM == 1 & 
                                HetVar_table_final$chrQTL_check==1) #10
eQTL_sumIM_check <- HetVar_table_final[idx_eQTL_sumIM_check,]
eQTL_sumIM_check_vars <- data.frame(Variants=paste(paste(eQTL_sumIM_check[,1],eQTL_sumIM_check[,2], sep=":"),
                                                   eQTL_sumIM_check[,3], sep="-"))
write.table(eQTL_sumIM_check_vars, "eQTL_sumIM_check_vars.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
eQTL_sumIM_check[,c(2,4,15)]

# expanding CpG IM to the entire peak (not just the summit)
sum(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$Peak_IM == 1 &
HetVar_table_final$chrQTL_check==1) #40
idx_eQTL_peakIM_check <- which(HetVar_table_final$Gene_eQTL==1 & HetVar_table_final$Peak_IM == 1 & 
                                HetVar_table_final$chrQTL_check==1) #40
eQTL_peakIM_check <- HetVar_table_final[idx_eQTL_peakIM_check,]
eQTL_peakIM_check_vars <- data.frame(Variants=paste(paste(eQTL_peakIM_check[,1],eQTL_peakIM_check[,2], sep=":"),
                                                    eQTL_peakIM_check[,3], sep="-"))
dim(eQTL_peakIM_check_vars)
eQTL_peakIM_check_vars_fil <- 
  data.frame(Variants=eQTL_peakIM_check_vars[!(eQTL_peakIM_check_vars$Variants %in%
                                                 eQTL_sumIM_check_vars$Variants),])
dim(eQTL_peakIM_check_vars_fil) # 31
write.table(eQTL_peakIM_check_vars_fil, "eQTL_peakIM_check_vars_fil.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
eQTL_peakIM_check[,c(2,4,15)]

# What about miRNAs?
sum(HetVar_table_final$miRNA_eQTL==1 & HetVar_table_final$Summit_IM == 1 &
      HetVar_table_final$chrQTL_check==1) #21
idx_miRNA_eQTL_sumIM_check <- which(HetVar_table_final$miRNA_eQTL==1 & HetVar_table_final$Summit_IM == 1 &
                                      HetVar_table_final$chrQTL_check==1) #21
length(idx_miRNA_eQTL_sumIM_check) #21
miRNA_eQTL_sumIM_check <- HetVar_table_final[idx_miRNA_eQTL_sumIM_check,]
miRNA_eQTL_sumIM_check_vars <- data.frame(Variants=paste(paste(miRNA_eQTL_sumIM_check[,1],miRNA_eQTL_sumIM_check[,2], sep=":"),
                                                         miRNA_eQTL_sumIM_check[,3], sep="-"))
write.table(miRNA_eQTL_sumIM_check_vars, "miRNA_eQTL_sumIM_check_vars.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
miRNA_eQTL_sumIM_check[,c(2,4,16)]


sum(HetVar_table_final$miRNA_eQTL==1 & HetVar_table_final$Peak_IM == 1 &
      HetVar_table_final$chrQTL_check==1) #79
idx_miRNA_eQTL_peakIM_check <- which(HetVar_table_final$miRNA_eQTL==1 & HetVar_table_final$Peak_IM == 1 &
                                      HetVar_table_final$chrQTL_check==1) #21
length(idx_miRNA_eQTL_peakIM_check) #79
miRNA_eQTL_peakIM_check <- HetVar_table_final[idx_miRNA_eQTL_peakIM_check,]
miRNA_eQTL_peakIM_check_vars <- data.frame(Variants=paste(paste(miRNA_eQTL_peakIM_check[,1],miRNA_eQTL_peakIM_check[,2], sep=":"),
                                                         miRNA_eQTL_peakIM_check[,3], sep="-"))

head(miRNA_eQTL_peakIM_check_vars)
dim(miRNA_eQTL_peakIM_check_vars) #79
miRNA_eQTL_peakIM_check_vars_fil <- 
  data.frame(Variants=miRNA_eQTL_peakIM_check_vars[!(miRNA_eQTL_peakIM_check_vars$Variants %in%
                                                 miRNA_eQTL_sumIM_check_vars$Variants),])
dim(miRNA_eQTL_peakIM_check_vars_fil) #58
write.table(miRNA_eQTL_peakIM_check_vars_fil, "miRNA_eQTL_peakIM_check_fil_vars.txt", append = F, quote = F,
            row.names = F, col.names = F, sep = "\t")
miRNA_eQTL_peakIM_check[,c(2,4,16)]


table(as.character(HetVar_table_final[HetVar_table_final$miRNA=="MIMAT0004608",7]))

colnames(HetVar_table_final)
head(HetVar_table_final)

### Make UpsetR plot 
library(UpSetR)

Gene_eQTLs<- readRDS("../../RNA-seq/Quant_2018/Genes_haplo_protlinc_sig2.rds")
head(Gene_eQTLs)
dim(Gene_eQTLs)

# miRNA 
miRNA_eQTL <- readRDS("../miRNA_QTLs/miRNA_haplo_eQTL_edit_matpat.rds")
head(miRNA_eQTL)
dim(miRNA_eQTL) # 47  25

chrQTL <- readRDS("../ATAC_summit_QTL_2/ATAC_haplo_chrQTL_edit_matpat.rds")
head(chrQTL)

chrQTL[chrQTL$V10=="Haplo_31",]


# merge them all on haplo type
Gene_eQTLs_haplo <- Gene_eQTLs[,1:2]
miRNA_eQTL_haplo <- miRNA_eQTL[,c(7,4)]
chrQTL_haplo <- chrQTL[,c(7,4)]
sum(table(chrQTL_haplo$V10)>0)

gene_mirna_haplo <- merge(Gene_eQTLs_haplo, miRNA_eQTL_haplo, by.x="V10", by.y="V10", all=TRUE)
head(gene_mirna_haplo)
dim(gene_mirna_haplo) # 363

gene_mirna_chr_haplo <- merge(gene_mirna_haplo, chrQTL_haplo, by.x="V10", by.y="V10", all=TRUE)
head(gene_mirna_chr_haplo)
dim(gene_mirna_chr_haplo) # 726   4
gene_mirna_chr_haplo <- gene_mirna_chr_haplo[!duplicated(gene_mirna_chr_haplo$V10),]
dim(gene_mirna_chr_haplo)

gene_mirna_chr_haplo2 <- gene_mirna_chr_haplo
head(gene_mirna_chr_haplo)
gene_mirna_chr_haplo2$caQTL <- 0
gene_mirna_chr_haplo2$caQTL[!is.na(gene_mirna_chr_haplo$V4)] <- 1
gene_mirna_chr_haplo2$gene_eQTL <- 0
gene_mirna_chr_haplo2$gene_eQTL[!is.na(gene_mirna_chr_haplo$V4.x)] <- 1
gene_mirna_chr_haplo2$miRNA_eQTL <- 0
gene_mirna_chr_haplo2$miRNA_eQTL[!is.na(gene_mirna_chr_haplo$V4.y)] <- 1
dim(gene_mirna_chr_haplo2)


gene_mirna_chr_haplo2 <- gene_mirna_chr_haplo2[,c(4:7)]
colnames(gene_mirna_chr_haplo2)[3:4] <- c("gene eQTL", "miRNA eQTL")

pdf(file = "Upset_haplotype_blocks_2018.pdf", family="ArialMT", width=6, height=4)
upset(gene_mirna_chr_haplo2, order.by="freq", empty.intersections="on",
      sets.x.label="Haplotype Blocks")
dev.off()
















