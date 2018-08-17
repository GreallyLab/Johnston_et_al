# Downsampling miRNA discovery
## Andrew D. Johnston 
### 03/05/18

We wanted to discern what depth of coverage was necessary for miRNA discovery. To do this we downsampled the uniquely mapped, non-duplicated, primary alignments of the merged sample sRNA datasets. The mapped reads are located in `~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/`. Additionally, all reads mapping to decoy and non-canonical chromosomes were excluded in the total read counts (chrEBV was included as a control).

Therefore, the total read count before exclusion was:

```bash
cd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/

module load samtools
samtools flagstat All.Aligned.sortedByCoord.out.filtered.primary.bam >  All.Aligned.sortedByCoord.out.filtered.primary.flagstat

more All.Aligned.sortedByCoord.out.filtered.primary.flagstat
 # 657,893,266 million reads are primary mapped. 

 # I need to piece together the genome after excluding the decoy, EBV chromosomes, and non-canonical chromosomes. 
mkdir Sim_discovery
cp split_primary_can/all.filtered.primary.chr16.bam Sim_discovery/
cp split_primary_can17/all.filtered.primary.chr17.bam Sim_discovery/

 # removing the EBV chr from the total number of reads
cd split_primary_can18up/

cp split_primary_can17/all.filtered.primary.chr18up.bam Sim_discovery/
mv All.Aligned.sortedByCoord.out.filtered.primary.REF_chrEBV.bam  ../

qsub -S /bin/bash -N merge_bam -l h_vmem=6G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.chr18up_noEBV.bam *REF*.bam
EOF 
mv ../All.Aligned.sortedByCoord.out.filtered.primary.REF_chrEBV.bam  .
mv all.filtered.primary.chr18up_noEBV.bam ../Sim_discovery

 # grabbing chr21
cd ../
cp split_primary_can21/All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.bam Sim_discovery/

 # finally create the final bam
cd Sim_discovery
qsub -S /bin/bash -N merge_bam -q highmem.q -l h_vmem=30G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.canchr.bam *.bam
samtools index all.filtered.primary.canchr.bam
samtools idxstats all.filtered.primary.canchr.bam > all.filtered.primary.canchr.idxstats.txt
samtools flagstat all.filtered.primary.canchr.bam > all.filtered.primary.canchr.flagstat.txt
EOF
```

There are 517,108,775 total reads. 

With this file, I plan to create various subsamples and run the miRDeep* software. 1 million reads, 2 mil, 5 mil, 10 mil, 20 mil, 50 mil, 100 mil, 250mil, 500 mil. 

Desired amount	Proportion
 1,000,000 	 0.0019338 
 2,000,000 	 0.0038677 
 5,000,000 	 0.0096691 
 10,000,000 	 0.0193383 
 25,000,000 	 0.0483457 
 50,000,000 	 0.0966915 
 100,000,000 	 0.1933829 
 250,000,000 	 0.4834573 
 500,000,000 	 0.9669146 

Making the subsamples:

```bash
cd Sim_discovery

qsub -S /bin/bash -N ss_bam_1 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss1
samtools view -s 0.0019338 -b all.filtered.primary.canchr.bam > ss1/all.filtered.primary.canchr.1mil.bam
EOF

qsub -S /bin/bash -N ss_bam_2 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss2
samtools view -s 0.0038677 -b all.filtered.primary.canchr.bam > ss2/all.filtered.primary.canchr.2mil.bam
EOF

qsub -S /bin/bash -N ss_bam_5 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss5
samtools view -s 0.0096691 -b all.filtered.primary.canchr.bam > ss5/all.filtered.primary.canchr.5mil.bam
EOF

qsub -S /bin/bash -N ss_bam_10 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss10
samtools view -s 0.0193383 -b all.filtered.primary.canchr.bam > ss10/all.filtered.primary.canchr.10mil.bam
EOF

qsub -S /bin/bash -N flagstat -l h_vmem=5G -cwd -j y << EOF
module load samtools
samtools flagstat all.filtered.primary.canchr.10mil.bam > all.filtered.primary.canchr.10mil.flagstat
EOF

qsub -S /bin/bash -N ss_bam_25 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss25
samtools view -s 0.0483457 -b all.filtered.primary.canchr.bam > ss25/all.filtered.primary.canchr.25mil.bam
EOF

qsub -S /bin/bash -N flagstat -l h_vmem=5G -cwd -j y << EOF
module load samtools
samtools flagstat all.filtered.primary.canchr.25mil.bam > all.filtered.primary.canchr.25mil.flagstat
EOF

qsub -S /bin/bash -N ss_bam_50 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss50
samtools view -s 0.0966915 -b all.filtered.primary.canchr.bam > ss50/all.filtered.primary.canchr.50mil.bam
EOF

qsub -S /bin/bash -N flagstat -l h_vmem=5G -cwd -j y << EOF
module load samtools
samtools flagstat all.filtered.primary.canchr.50mil.bam > all.filtered.primary.canchr.50mil.flagstat
EOF

qsub -S /bin/bash -N ss_bam_100 -l h_vmem=20G -cwd -j y << EOF
module load samtools
mkdir ss100
samtools view -s 0.1933829 -b all.filtered.primary.canchr.bam > ss100/all.filtered.primary.canchr.100mil.bam
EOF

qsub -S /bin/bash -N flagstat -l h_vmem=5G -cwd -j y << EOF
module load samtools
samtools flagstat all.filtered.primary.canchr.100mil.bam > all.filtered.primary.canchr.100mil.flagstat
EOF

qsub -S /bin/bash -N ss_bam_250 -l h_vmem=30G -cwd -j y << EOF
module load samtools
mkdir ss250
samtools view -s 0.4834573 -b all.filtered.primary.canchr.bam > ss250/all.filtered.primary.canchr.250mil.bam
EOF


qsub -S /bin/bash -N flagstat -l h_vmem=5G -cwd -j y << EOF
module load samtools
samtools flagstat all.filtered.primary.canchr.250mil.bam > all.filtered.primary.canchr.250mil.flagstat
EOF

qsub -S /bin/bash -N ss_bam_500 -l h_vmem=30G -cwd -j y << EOF
module load samtools
mkdir ss500
samtools view -s 0.9669146 -b all.filtered.primary.canchr.bam > ss500/all.filtered.primary.canchr.500mil.bam
EOF

qsub -S /bin/bash -N flagstat -l h_vmem=15G -cwd -j y << EOF
module load samtools
samtools flagstat all.filtered.primary.canchr.500mil.bam > all.filtered.primary.canchr.500mil.flagstat
EOF

```

I won't have to split the smaller subsamples, but may have to split the chromosomes up for 250mil and 500mil. Let's see. I am running the command line version of MrDeep*, which is a Java program that requires a ton of memory. 

```bash
cd ~/Programs/MDS_command_line_v37/MDS_command_line

qsub -S /bin/bash -N ss1_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss1/all.filtered.primary.canchr.1mil.bam
EOF

qsub -S /bin/bash -N ss2_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss2/all.filtered.primary.canchr.2mil.bam
EOF

qsub -S /bin/bash -N ss5_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss5/all.filtered.primary.canchr.5mil.bam
EOF

qsub -S /bin/bash -N ss10_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss10/all.filtered.primary.canchr.10mil.bam
EOF


qsub -S /bin/bash -N ss25_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss25/all.filtered.primary.canchr.25mil.bam
EOF

qsub -S /bin/bash -N ss50_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss50/all.filtered.primary.canchr.50mil.bam
EOF

qsub -S /bin/bash -N ss100_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss100/all.filtered.primary.canchr.100mil.bam
EOF

qsub -S /bin/bash -N ss250_MD -q highmem.q -l h_vmem=1305G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss250/all.filtered.primary.canchr.250mil.bam
EOF
 
 # Our HPC does not have the memory capacity to handle 500 million reads with miRDeep*
qsub -S /bin/bash -N ss500_MD -q highmem.q -l h_vmem=1306G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss500/all.filtered.primary.canchr.500mil.bam
EOF

 # it got up to chr 17; awk '{print $3}' all.filtered.primary.canchr.500mil.result | sort | uniq -c 
 # therefore I will split up the reads from 1-16 and 17 and up
cd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss500

qsub -S /bin/bash -N split_ss500 -l h_vmem=30G -cwd -j y << EOF
module load bamtools
bamtools split -in all.filtered.primary.canchr.500mil.bam -reference
EOF

mkdir 16down
mkdir 17up

 #placed the respective chromosomes in each folder
cd 16down
qsub -S /bin/bash -N merge_16 -l h_vmem=30G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.chr16down.bam *.bam
EOF

cd ../17up
qsub -S /bin/bash -N merge_17 -l h_vmem=30G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.chr17up.bam *.bam
EOF

cd ~/Programs/MDS_command_line_v37/MDS_command_line

qsub -S /bin/bash -N ss500_16_MD -q highmem.q -l h_vmem=1306G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss500/16down/all.filtered.primary.chr16down.bam
EOF

qsub -S /bin/bash -N ss500_17_MD -q highmem.q -l h_vmem=1306G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss500/17up/all.filtered.primary.chr17up.bam
EOF

 # 17 up didn't work. So we'll have to perform similar splitting as done previously
cd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss500/17up/
 # need to remove MIR155
module load samtools
samtools index all.filtered.primary.canchr.500mil.REF_chr21.bam
samtools view all.filtered.primary.canchr.500mil.REF_chr21.bam chr21:25573980-25574044 | wc -l # 3,441,3258 reads (which is the issue)
samtools view -b all.filtered.primary.canchr.500mil.REF_chr21.bam chr21:1-25573980 chr21:25574044-46709983 > all.filtered.primary.canchr.500mil.REF_chr21minusMIR155.bam

 # need to remove MIR21 as well

qsub -S /bin/bash -N spl_chr17 -q highmem.q -l h_vmem=30G -cwd -j y << EOF
module load samtools
samtools index all.filtered.primary.canchr.500mil.REF_chr17.bam
samtools view all.filtered.primary.canchr.500mil.REF_chr17.bam chr17:59841266-59841337 | wc -l
 # 31,303,353 million which is microRNA 21.. which is also breaking the program. 
samtools view -b all.filtered.primary.canchr.500mil.REF_chr17.bam chr17:1-59841266 chr17:59841337-83257441 > all.filtered.primary.canchr.500mil.REF_chr17minusMIR21.bam
EOF

 # moved the unused bams to unused_bams/

qsub -S /bin/bash -N merge_17 -l h_vmem=30G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.chr17upMinus15521.bam *.bam
EOF

cd ~/Programs/MDS_command_line_v37/MDS_command_line

qsub -S /bin/bash -N ss500_17.2_MD -q highmem.q -l h_vmem=1306G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1200g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/ss500/17up/all.filtered.primary.chr17upMinus15521.bam
EOF
 # worked

```

After calling the miRDeep on different read depth, now we need to filter the novel miRNAs in a smart way... but we really only care if the novel miRNAs called by the downsample are the same as our curated one (from the max amount of depth). 

```bash
 # get the total number of novel miRNAs found (before filtering)
cd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Sim_discovery/
module load bedtools2

 # ss1
cd ss1
grep novel all.filtered.primary.canchr.1mil.result | wc -l # 38
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.1mil.result | wc -l # 3
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.1mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_1mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.1mil.result | paste - AllStart_stop_canchr_1mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_1mil.bed

bedtools merge -i All_novel_canchr_1mil.bed | wc -l #35
bedtools merge -i All_novel_canchr_1mil.bed | bedtools intersect -a ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #9
bedtools merge -i All_novel_canchr_1mil.bed | bedtools intersect -a ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #17

awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.1mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_1mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.1mil.result | paste - FilStart_stop_canchr_1mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_1mil.bed

bedtools merge -i Fil_novel_canchr_1mil.bed | wc -l #3
bedtools merge -i Fil_novel_canchr_1mil.bed | bedtools intersect -a ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #2
bedtools merge -i Fil_novel_canchr_1mil.bed | bedtools intersect -a ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #3
```

```bash
 # ss2
cd ../ss2
grep novel all.filtered.primary.canchr.2mil.result | wc -l # 70
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.2mil.result | wc -l # 4
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.2mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_2mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.2mil.result | paste - AllStart_stop_canchr_2mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_2mil.bed

bedtools merge -i All_novel_canchr_2mil.bed | wc -l #67
bedtools merge -i All_novel_canchr_2mil.bed | bedtools intersect -a ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #11
bedtools merge -i All_novel_canchr_2mil.bed | bedtools intersect -a ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #24
bedtools intersect -a All_novel_canchr_2mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0

awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.2mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_2mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.2mil.result | paste - FilStart_stop_canchr_2mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_2mil.bed

bedtools merge -i Fil_novel_canchr_2mil.bed | wc -l #4
bedtools merge -i Fil_novel_canchr_2mil.bed | bedtools intersect -a ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #3
bedtools merge -i Fil_novel_canchr_2mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #4
bedtools intersect -a Fil_novel_canchr_2mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```

```bash
 # ss5
cd ../ss5
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.5mil.result | wc -l # 109
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.5mil.result | wc -l # 5

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.5mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_5mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.5mil.result | paste - AllStart_stop_canchr_5mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_5mil.bed

bedtools merge -i All_novel_canchr_5mil.bed | wc -l #104
bedtools merge -i All_novel_canchr_5mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #14
bedtools merge -i All_novel_canchr_5mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #30
bedtools intersect -a All_novel_canchr_5mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.5mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_5mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.5mil.result | paste - FilStart_stop_canchr_5mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_5mil.bed

bedtools merge -i Fil_novel_canchr_5mil.bed | wc -l #5
bedtools intersect -a Fil_novel_canchr_5mil.bed -b ../Pooled_novel_miRNA_merge.bed -u -wb | wc -l #4
bedtools intersect -a Fil_novel_canchr_5mil.bed -b ../Combined_novel_miRNA_merge.bed -u -wb | wc -l #5
bedtools intersect -a Fil_novel_canchr_5mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```


```bash
 # ss10
cd ../ss10
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.10mil.result | wc -l # 163
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.10mil.result | wc -l # 4

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.10mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_10mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.10mil.result | paste - AllStart_stop_canchr_10mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_10mil.bed

bedtools merge -i All_novel_canchr_10mil.bed | wc -l #155
bedtools merge -i All_novel_canchr_10mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #21
bedtools merge -i All_novel_canchr_10mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #38
bedtools intersect -a All_novel_canchr_10mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.10mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_10mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.10mil.result | paste - FilStart_stop_canchr_10mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_10mil.bed

bedtools merge -i Fil_novel_canchr_10mil.bed | wc -l #4
bedtools intersect -a Fil_novel_canchr_10mil.bed -b ../Pooled_novel_miRNA_merge.bed -u -wb | wc -l #3
bedtools intersect -a Fil_novel_canchr_10mil.bed -b ../Combined_novel_miRNA_merge.bed -u -wb | wc -l #4
bedtools intersect -a Fil_novel_canchr_10mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```

```bash
 # ss25
cd ../ss25
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.25mil.result | wc -l # 325
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.25mil.result | wc -l # 6

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.25mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_25mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.25mil.result | paste - AllStart_stop_canchr_25mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_25mil.bed

bedtools merge -i All_novel_canchr_25mil.bed | wc -l #307
bedtools merge -i All_novel_canchr_25mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #42
bedtools merge -i All_novel_canchr_25mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #61
bedtools intersect -a All_novel_canchr_25mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.25mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_25mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.25mil.result | paste - FilStart_stop_canchr_25mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_25mil.bed

bedtools merge -i Fil_novel_canchr_25mil.bed | wc -l #6
bedtools intersect -a Fil_novel_canchr_25mil.bed -b ../Pooled_novel_miRNA_merge.bed -u -wb | wc -l #4
bedtools intersect -a Fil_novel_canchr_25mil.bed -b ../Combined_novel_miRNA_merge.bed -u -wb | wc -l #6
bedtools intersect -a Fil_novel_canchr_25mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```

```bash
 # ss50
cd ../ss50
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.50mil.result | wc -l # 578
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.50mil.result | wc -l # 13

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.50mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_50mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.50mil.result | paste - AllStart_stop_canchr_50mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_50mil.bed 

bedtools merge -i All_novel_canchr_50mil.bed | wc -l #548
bedtools merge -i All_novel_canchr_50mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #77
bedtools merge -i All_novel_canchr_50mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #96
bedtools intersect -a All_novel_canchr_50mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #1

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.50mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_50mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.50mil.result | paste - FilStart_stop_canchr_50mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_50mil.bed

bedtools merge -i Fil_novel_canchr_50mil.bed | wc -l #548
bedtools intersect -a Fil_novel_canchr_50mil.bed -b ../Pooled_novel_miRNA_merge.bed -u -wb | wc -l #11
bedtools intersect -a Fil_novel_canchr_50mil.bed -b ../Combined_novel_miRNA_merge.bed -u -wb | wc -l #13
bedtools intersect -a Fil_novel_canchr_50mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```

```bash
 # ss100
cd ../ss100
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.100mil.result | wc -l # 978
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.100mil.result | wc -l # 28

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.100mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_100mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.100mil.result | paste - AllStart_stop_canchr_100mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_100mil.bed 

bedtools merge -i All_novel_canchr_100mil.bed | wc -l #922
bedtools merge -i All_novel_canchr_100mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #140
bedtools merge -i All_novel_canchr_100mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u -wb | wc -l #160
bedtools intersect -a All_novel_canchr_100mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #1

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.100mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_100mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.100mil.result | paste - FilStart_stop_canchr_100mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_100mil.bed

bedtools merge -i Fil_novel_canchr_100mil.bed | wc -l #28
bedtools intersect -a Fil_novel_canchr_100mil.bed -b ../Pooled_novel_miRNA_merge.bed -u -wb | wc -l #22
bedtools intersect -a Fil_novel_canchr_100mil.bed -b ../Combined_novel_miRNA_merge.bed -u -wb | wc -l #25
bedtools intersect -a Fil_novel_canchr_100mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```

```bash
 # ss250
cd ../ss250

awk '{print $3}' all.filtered.primary.canchr.250mil.result | sort | uniq -c
 # all of the chromosomes are represented, so it's all fine
 
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.250mil.result | wc -l # 2186
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.250mil.result | wc -l # 76

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.250mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_250mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.250mil.result | paste - AllStart_stop_canchr_250mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_250mil.bed 

bedtools merge -i All_novel_canchr_250mil.bed | wc -l #2051
bedtools merge -i All_novel_canchr_250mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l #	
bedtools merge -i All_novel_canchr_250mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l #174
bedtools intersect -a All_novel_canchr_250mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #1

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.250mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_250mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.250mil.result | paste - FilStart_stop_canchr_250mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_250mil.bed

bedtools merge -i Fil_novel_canchr_250mil.bed | wc -l #75
bedtools merge -i Fil_novel_canchr_250mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u -wb | wc -l #64
bedtools merge -i Fil_novel_canchr_250mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u -wb | wc -l #66
bedtools intersect -a Fil_novel_canchr_250mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #0
```

For 500 million reads, I need to combine the 16down and 17up results files

```bash
 # ss500
cd ../ss500

awk 'NR>1' 17up/all.filtered.primary.chr17upMinus15521.result | cat 16down/all.filtered.primary.chr16down.result - > all.filtered.primary.canchr.500mil.result

awk '{print $3}' all.filtered.primary.canchr.500mil.result | sort | uniq -c
 # all of the chromosomes are represented, so it's all fine
 
 # how many novel miRNAs are called?
grep novel all.filtered.primary.canchr.500mil.result | wc -l # 2186
 # novel after using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.canchr.500mil.result | wc -l # 76

 # creating a bed of all novel miRNAs
awk '$1~"novel" {print $7}' all.filtered.primary.canchr.500mil.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_canchr_500mil.txt
awk '$1~"novel" {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.500mil.result | paste - AllStart_stop_canchr_500mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > All_novel_canchr_500mil.bed 

bedtools merge -i All_novel_canchr_500mil.bed | wc -l #3774
bedtools merge -i All_novel_canchr_500mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u | wc -l # 163
bedtools merge -i All_novel_canchr_500mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u | wc -l # 180
bedtools intersect -a All_novel_canchr_500mil.bed -b ../novelMiR_550.bed -u -wb | wc -l #1

 # creating a bed using same filters as the pooled analysis
awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.canchr.500mil.result | awk -F "-" '{print $1"\t"$2}' > FilStart_stop_canchr_500mil.txt
awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.canchr.500mil.result | paste - FilStart_stop_canchr_500mil.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n > Fil_novel_canchr_500mil.bed

bedtools merge -i Fil_novel_canchr_500mil.bed | wc -l #162
bedtools merge -i Fil_novel_canchr_500mil.bed | bedtools intersect -a  ../Pooled_novel_miRNA_merge.bed -b - -u -wb | wc -l #158
bedtools merge -i Fil_novel_canchr_500mil.bed | bedtools intersect -a  ../Combined_novel_miRNA_merge.bed -b - -u -wb | wc -l #158
bedtools intersect -a Fil_novel_canchr_500mil.bed -b ../novelMiR_550.bed -u -wb | wc -l # 1
```
 I need to develop a list of miRNAs found that are not in Gencode and HSA (but in fantom and noncode)


Because the pooled analysis and the individual analysis included all chromosomes and not just canonical, I need to figure out how many novel miRNA are in just canon chr.

```bash
wc -l Combined_novel_miRNA.bed 
 # 221 Combined_novel_miRNA.bed
 
sort -k 1,1 -k2,2n Combined_novel_miRNA.bed  | bedtools merge -i - | wc -l # 220
sort -k 1,1 -k2,2n Combined_novel_miRNA.bed  | bedtools merge -i - > Combined_novel_miRNA_merge.bed

awk '{print $1}' Combined_novel_miRNA_merge.bed | sort | uniq -c 
```
These are the only non canon chrs with novel miRNAS 
 	  6 chrEBV
      1 chrJTFH01000007.1
      1 chrJTFH01000009.1

Therefore, there's 212 canon chr 

```bash
wc -l Pooled_novel_miRNA.bed 
 # 170 Pooled_novel_miRNA.bed
 
sort -k 1,1 -k2,2n Pooled_novel_miRNA.bed  | bedtools merge -i - | wc -l # 168
sort -k 1,1 -k2,2n Pooled_novel_miRNA.bed  | bedtools merge -i - > Pooled_novel_miRNA_merge.bed

awk '{print $1}' Pooled_novel_miRNA_merge.bed | sort | uniq -c 
```
These are the only non canon chrs with novel miRNAS 
 	  3 chrEBV
Therefore, there's 165 canon chr in pooled analysis

