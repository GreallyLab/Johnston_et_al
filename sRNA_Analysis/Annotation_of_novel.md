After searching for novel miRNAs by passing all sRNAseq reads into miRDeep*, I needed to filter for the most likely miRNAs as well as annotate them.

I filtered the individuals blocks of the genome for which novel miRNAs were sought. 

working in `/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can`

```bash
 # 16 and down
 # number of novel?
grep novel all.filtered.primary.chr16.result | wc -l 
 # 2748

 # number of known
grep hsa all.filtered.primary.chr16.result | wc -l 
 # 463

 # median score of known and novel
awk '$1~"hsa"{print $2}' all.filtered.primary.chr16.result | sort -n | awk -f median.awk 
 # 2.12 (3.785) (15.33)
awk '$1~"novel" {print $2}' all.filtered.primary.chr16.result | sort -n | awk -f median.awk 
 # -1.675 (-1.60) (-2.08)

 # how many novel pass the filter?
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.chr16.result | wc -l #104 (59) novel miRNA

 awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.chr16.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr16.txt

 awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.chr16.result | paste - AllStart_stop_chr16.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr16.bed # 59

 cp All_novel_chr16.bed ../Novel_miRNA/
```


 ######################### 17 1half #############################
working in `~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can17`
```bash

cp ../split_primary_can/median.awk .

 # number of novel?
grep novel all.filtered.primary.chr17.1half.result | wc -l 
 # 138

 # number of known
grep hsa all.filtered.primary.chr17.1half.result | wc -l 
 # 26

 # median score of known and novel
awk '$1~"hsa"{print $2}' all.filtered.primary.chr17.1half.result | sort -n | awk -f median.awk 
 # 25.965
awk '$1~"novel" {print $2}' all.filtered.primary.chr17.1half.result | sort -n | awk -f median.awk 
 # -0.47

 # how many novel pass the filter?
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.chr17.1half.result | wc -l #8 

awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.chr17.1half.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr17_1half.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.chr17.1half.result | paste - AllStart_stop_chr17_1half.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr17_1half.bed # 8

cp All_novel_chr17_1half.bed ../Novel_miRNA/


 ##################### 17 2half ######################
 # number of novel?
grep novel all.filtered.primary.chr17.2half.minusMIR21.result| wc -l 
 # 192

 # number of known
grep hsa all.filtered.primary.chr17.2half.minusMIR21.result| wc -l 
 # 25

 # median score of known and novel
awk '$1~"hsa"{print $2}' all.filtered.primary.chr17.2half.minusMIR21.result| sort -n | awk -f median.awk 
 # 25.965
awk '$1~"novel" {print $2}' all.filtered.primary.chr17.2half.minusMIR21.result| sort -n | awk -f median.awk 
 # -1.75

 # how many novel pass the filter?
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.chr17.2half.minusMIR21.result| wc -l #8 

awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.chr17.2half.minusMIR21.result| awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr17_2half.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.chr17.2half.minusMIR21.result| paste - AllStart_stop_chr17_2half.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr17_2half.bed # 8

cp All_novel_chr17_2half.bed ../Novel_miRNA/
```

Then for 18 and up
 ######################### 18andUp ###########################
working in  `pwd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can18up`
```bash
cp ../split_primary_can/median.awk .

 # number of novel?
grep novel all.filtered.primary.chr18up.result | wc -l 
 # 1050

 # number of known
grep hsa all.filtered.primary.chr18up.result| wc -l 
 # 184

 # median score of known and novel
awk '$1~"hsa"{print $2}' all.filtered.primary.chr18up.result| sort -n | awk -f median.awk 
 # 2.08
awk '$1~"novel" {print $2}' all.filtered.primary.chr18up.result| sort -n | awk -f median.awk 
 # -1.71

 # how many novel pass the filter?
awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.chr18up.result| wc -l # 49

awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.chr18up.result| awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr18Up.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.chr18up.result| paste - AllStart_stop_chr18Up.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr18up.bed # 8

cp All_novel_chr18up.bed ../Novel_miRNA/
 ########################## 21 1st half ####################### 
 # pwd /17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can21
cp ../split_primary_can/median.awk .

 # number of novel?
grep novel All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result | wc -l 
 # 12

 # number of known
grep hsa All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result| wc -l 
 # 4

 # median score of known and novel
awk '$1~"hsa"{print $2}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result| sort -n | awk -f median.awk 
 # 4518.64
awk '$1~"novel" {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result| sort -n | awk -f median.awk 
 # -1.155

 # how many novel pass the filter?
awk '$1~"novel" && $6>50 && $2>10 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result | wc -l # 1

awk '$1~"novel" && $6>50 && $2>10 {print $7}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result| awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr21_1half.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.result| paste - AllStart_stop_chr21_1half.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr21_1half.bed # 8

cp All_novel_chr21_1half.bed ../Novel_miRNA/

 ############################## 21 2half#######################

 # number of novel?
grep novel All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result | wc -l 
 # 44

 # number of known
grep hsa All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result| wc -l 
 # 2

 # median score of known and novel
awk '$1~"hsa"{print $2}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result| sort -n | awk -f median.awk 
 # 281.565
awk '$1~"novel" {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result| sort -n | awk -f median.awk 
 # -1.22

 # how many novel pass the filter?
awk '$1~"novel" && $6>50 && $2>10 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result| wc -l # 1

awk '$1~"novel" && $6>50 && $2>10 {print $7}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result| awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr21_2half.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.result| paste - AllStart_stop_chr21_2half.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr21_2half.bed # 8

cp All_novel_chr21_2half.bed ../Novel_miRNA/
```


Now I combined all of the novel miRNAs discovered by pooling all of the sRNAseq libraries. 

```bash
 # combining all of the novel miRNAs from the pooled reads 
cd ../Novel_miRNA/

cat *.bed > Pooled_novel_miRNA.bed
wc -l Pooled_novel_miRNA.bed
 # 170 Pooled_novel_miRNA.bed

 # I will then check how many novel miRNAs discovered using individual samples were discovered by the pooling technique
 # compared to the discovery when miRNAs are discovered within each individual only. 
module load bedtools2
 # edited ../Merged_novel_fillist.bed to a bed file with only miRNAs passing the score of 10 threshold, creating a new file Individual_novel_miRNA.bed

bedtools intersect -a Individual_novel_miRNA.bed -b Pooled_novel_miRNA.bed -wo | wc -l
 # 27 miRNAs were found by both techniques.

bedtools intersect -a Individual_novel_miRNA.bed -b Pooled_novel_miRNA.bed -wo > Pooled_individual_overlap.txt

 # Then I will combine the individual and pooled technique novel miRNAs see Combine_novel.R
wc -l Combined_novel_miRNA.bed
 # 221 Combined_novel_miRNA.bed
```

I want to incorporate the novel miRNAs into the quantification via featureCounts, so I made a gtf file

```bash
 ###### Generate counts for all of the novel miRNAs ##########

 # I need to make a GTF file for featureCounts to use. 

 # had to manually take out the chr in front of the non-chr contigs such as GL000220.1 to get the fasta seqeunce (because index fasta didn't have them).
cp Combined_novel_miRNA.bed Combined_novel_miRNA_chrfix.bed

less miRNA_hg38_EBV_decoy_sort.gtf

awk '{print $4}' Combined_novel_miRNA_chrfix.bed | awk '{sub($1, "\"&\""); print}' > Combined_novel_miRNA_chrfix_names.txt

awk '{print $1"\t.\tmiRNA\t"$2"\t"$3"\t"$5"\t"$6"\t.\tID"}' Combined_novel_miRNA_chrfix.bed | paste -d " " - Combined_novel_miRNA_chrfix_names.txt > Combined_novel_miRNA_chrfix.gtf

 #wc -l Combined_novel_miRNA_chrfix.gtf
221 Combined_novel_miRNA_chrfix.gtf

cat /home/greally-lab/indexes/hg38/miRNA/hsa.GTF Combined_novel_miRNA_chrfix.gtf | sort -k 1,1 -k 4,2n > miRNA_hg38_EBV_decoy_plus_novel.gtf

mkdir Count_miRNA_novel
 ## count the miRNAs using featureCounts - not using primary because I want the most stringent counting for the miRNAs. I will use the unfiltered bam file (not filtered for 5' soft clipping) because I'm more curious about counts not the exact positioning of the miRNA transcript

for f1 in ../*Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -d '/' -f2 | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N ${SAMPLE}_featureCounts -j y -cwd -q highmem.q -pe smp-highmem 2 -l h_vmem=5.6G << EOF
module load subread/1.5.0-p1/gcc.4.4.7 
featureCounts -T 2 -F GTF -t miRNA -g ID -a miRNA_hg38_EBV_decoy_plus_novel.gtf -o Count_miRNA_novel/${SAMPLE}_miRNA.txt ${f1}
EOF
done
```

I wanted to assess the novelty of the discovered miRNAs; therefore I checked several sRNA/miRNA databases to check whether or not my miRNA were found in these databases. I also wanted to annotate the known miRNAs (miRBase) as well. 

```bash
 # make a bed of primary and mature 
awk '{print $3}' /home/greally-lab/indexes/hg38/miRNA/hsa.GTF | sort | uniq -c 
 # 2813 miRNA
 # 1881 miRNA_primary_transcript
awk '$3 == "miRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' /home/greally-lab/indexes/hg38/miRNA/hsa.GTF | tr -d '";' > known_miRNA.bed
wc -l known_miRNA.bed
 #2813
 # ok so let's combine the bed files. 
cat known_miRNA.bed Combined_novel_miRNA_chrfix.bed > All_miRNA_chrfix.bed
 # Now let's go through all of the filtering steps.
wc -l All_miRNA_chrfix.bed
 # 3,034 All_miRNA_chrfix.bed
```


1. Rfam database v.12

I had to align my miRNAs to the Rfam database.

```bash
 # Now I'll see which are within Rfam
 #downloaded the v.12 Rfam fasta sequence to /home/greally-lab/indexes/
cd /home/greally-lab/indexes/

 # Build reference genome database
module load gmap/2016-09-23/gcc.4.4.7 
qsub -S /bin/bash -N Build_Rfam_gmap -j y -cwd -l h_vmem=40G << EOF
module load gmap/2016-09-23/gcc.4.4.7 
gmap_build -D Rfam_assembly -g -k 15 -d Rfam RFam_all.fa.gz
EOF

cd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Novel_miRNA
cp Combined_novel_miRNA.bed Combined_novel_miRNA_chrfix.bed

 # had to manually take out the chr in front of the non-chr contigs such as GL000220.1 to get the fasta seqeunce (because index fasta didn't have them).

module load bedtools2 
bedtools getfasta -name -fi /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa -bed All_miRNA_chrfix.bed -fo All_miRNA_chrfix.fasta

 # Alignment
qsub -S /bin/bash -N All_RFam_gmap -j y -cwd -pe smp 4 -l h_vmem=5G << EOF
module load gmap/2016-09-23/gcc.4.4.7 
module load samtools
gmap -D /home/greally-lab/indexes/Rfam_dir/Rfam_assembly/Rfam  --no-chimeras  --nosplicing -d Rfam -f samse -t 4 All_miRNA_chrfix.fasta | samtools view -Shb - | samtools sort - All_RFam_alignment
samtools index All_RFam_alignment.bam
samtools flagstat All_RFam_alignment.bam	
 # to see which families are involed, this needs to be filtered. 
samtools idxstats All_RFam_alignment.bam > All_RFam_alignment.idxstats.txt	
EOF

samtools flagstat All_RFam_alignment.bam
 # 4102 + 0 in total (QC-passed reads + QC-failed reads)
 # 1406 + 0 mapped (34.28%:-nan%)

 # 3034 sequences became 4102 because of some miRNAs aligning to multiple RFams, see below. 3034+1068 secondary alignments

samtools view All_RFam_alignment.bam | awk '$3!="*"' | awk '{print $1}' | sort | uniq -c | wc -l # 3034 seqeunces, 338 sequences overlapped with Rfam - that's 11.1%

samtools view All_RFam_alignment.bam | awk '$3!="*"' | awk '{print $1}' | sort | uniq -c > All_inter_RFam.txt 
```

2. filter by GenCode sRNA and miRNA

```bash
bedtools intersect -a All_miRNA_chrfix.bed -b ../gencode.v24.primary_assembly.annotation.sRNA.bed -wo > All_miRNA_inter_gencode_sRNA.txt
wc -l All_miRNA_inter_gencode_sRNA.txt # 227
awk '{print $11}' All_miRNA_inter_gencode_sRNA.txt | sort | uniq -c
 # 	 174 antisense
 #      4 rRNA
 #      4 scaRNA
 #     45 snoRNA

bedtools intersect -a All_miRNA_chrfix.bed -b gencode.v24.primary_assembly.annotation.miRNA.bed -wo | wc -l # 2734

bedtools intersect -a All_miRNA_chrfix.bed -b gencode.v24.primary_assembly.annotation.miRNA.bed -wo >  All_miRNA_inter_Gencode_miRNA.txt

 # Filter by Gencode tRNA
awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$6"\t"$8}' rRNAandtRNA_GenCode.gtf | tr -d '";' > gencode.v24.primary_assembly.annotation.tRNAnrRNA.bed

 cp /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/gencode.v24.primary_assembly.annotation.tRNAnrRNA.bed .
awk '$3~"gene" {print $12}' gencode.v24.primary_assembly.annotation.EBV.gtf | sort | uniq -c

bedtools intersect -a All_miRNA_chrfix.bed -b /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/gencode.v24.primary_assembly.annotation.tRNAnrRNA.bed -wo > All_miRNA_inter_gencode_rRNAntRNA.txt
wc -l All_miRNA_inter_gencode_rRNAntRNA.txt # 20
```

3. filter by fRNAdb 

```bash
 # downloaded all of the sequences in the fRNAdb  from http://togodb.biosciencedbc.jp/togodb/view/frnadb_summary#en
 # fRNAdb.csv in /17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/Novel_miRNA
wc -l fRNAdb.csv
 # 126985
awk 'BEGIN {FS=","} $5 ~ "Homo" {print}' fRNAdb.csv > fRNAdb_human.csv 
 # 124549

 # make a fasta file that I can gmap align to.
awk 'BEGIN {FS=","} {print ">"$1"-"$4"-"$7"\n"$8}' fRNAdb_human.csv | fold -w 80 > fRNAdb_human.fasta

 # Build reference genome database
module load gmap/2016-09-23/gcc.4.4.7 
qsub -S /bin/bash -N Build_fRNAdb_gmap -j y -cwd -q highmem.q -l h_vmem=50G << EOF
module load gmap/2016-09-23/gcc.4.4.7 
gmap_build -D fRNAdb_gmap -k 15 -d fRNAdb fRNAdb_human.fasta
EOF

 # Alignment
qsub -S /bin/bash -N All_fRNA_gmap -j y -cwd -pe smp 4 -l h_vmem=5G << EOF
module load gmap/2016-09-23/gcc.4.4.7 
module load samtools
gmap -D fRNAdb_gmap/fRNAdb --no-chimeras --nosplicing -d fRNAdb -f samse -t 4 All_miRNA_chrfix.fasta | samtools view -Shb - | samtools sort - All_fRNA_alignment
samtools index All_fRNA_alignment.bam
samtools flagstat All_fRNA_alignment.bam
 # 4325 + 0 in total (QC-passed reads + QC-failed reads)
 # 2206 + 0 mapped (51.01%:-nan%)
samtools idxstats All_fRNA_alignment.bam > All_fRNA_alignment.idxstats
EOF

 # 3034 seqeunces, 1291 repeats, 915 sequences overlapped with fRNAdb - that's 30.6%

samtools view All_fRNA_alignment.bam | awk '$3!="*"' | awk '{print $1}' | sort | uniq -c | wc -l # 915

samtools view All_fRNA_alignment.bam | awk '$3!="*"' | awk '{print $1}' | sort | uniq -c > All_miRNA_inter_fRNA.txt 

samtools view All_fRNA_alignment.bam | awk '$3!="*"' | awk '{print $3}' | sort | uniq -c > All_miRNA_inter_fRNA_annot.txt 

samtools view All_fRNA_alignment.bam | awk '$3!="*"' | awk '{print $1"\t"$3}' > All_miRNA_inter_fRNA_and_annot.txt 
```

4. filter by noncode 2016

```bash
 # http://www.noncode.org/download.php
 # downloaded NONCODE2016_human.fa.gz

 # Build reference genome database
module load gmap/2016-09-23/gcc.4.4.7 
qsub -S /bin/bash -N Build_Noncode_gmap -j y -cwd -q highmem.q -l h_vmem=50G << EOF
module load gmap/2016-09-23/gcc.4.4.7 
gmap_build -D Noncode_gmap -k 15 -d Noncode NONCODE2016_human.fa
EOF

 # Alignment
qsub -S /bin/bash -N All_Noncode_gmap -j y -cwd -pe smp 4 -l h_vmem=5G << EOF
module load gmap/2016-09-23/gcc.4.4.7 
module load samtools
gmap -D Noncode_gmap/Noncode --no-chimeras --nosplicing -d Noncode -f samse -t 4 All_miRNA_chrfix.fasta | samtools view -Shb - | samtools sort - All_Noncode_alignment
samtools index All_Noncode_alignment.bam
samtools flagstat All_Noncode_alignment.bam	
 # 3589 + 0 in total (QC-passed reads + QC-failed reads)
 # 994 + 0 mapped (27.70%:-nan%)
samtools idxstats All_Noncode_alignment.bam > All_Noncode_alignment.idxstats
EOF

 # 3034 seqeunces, 555 repeats,therefore 439 different sequences overlapped with noncode - that's 14.4%

samtools view All_Noncode_alignment.bam | awk '$3!="*"' | awk '{print $1}' | sort | uniq -c | wc -l 
 # 439

samtools view All_Noncode_alignment.bam | awk '$3!="*"' | awk '{print $1}' | sort | uniq -c > All_miRNA_inter_noncode.txt 
```

5. filter out blacklisted regions and poor mappability 


```bash 
 # I should also get rid of blacklisted sites 
 # I downloaded blacklisted hg38 sites from anshul Kundaje 
 # http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/

module load bedtools2
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/hg38.blacklist.bed -wo | less
 # none are in the blacklisted regions
```

6. filter out miRNA in FANTOM

```bash
 # downloaded the new fantom data from 
 # into ~/Genome_annotations_hg38/
awk '$3=="gene"' FANTOM_CAT.lv3_robust.gtf | wc -l  
 # 59,110 genes

awk '$3=="gene" && $12~"short_RNA"' FANTOM_CAT.lv3_robust.gtf | wc -l  
 # 2,678 short RNA genes

awk '$3=="gene" && $12~"short_RNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$28"\t.\t"$16}' FANTOM_CAT.lv3_robust.gtf | tr -d '";' > FANTOM_lv3_shortRNA.bed  

awk '$3=="gene" && $12~"short_RNA" {print $16}' FANTOM_CAT.lv3_robust.gtf | sort | uniq -c
 #	  915 "longer_than_100nt";
 #     594 "miRNA";
 #       2 "Mt_rRNA";
 #      22 "Mt_tRNA";
 #      40 "rRNA";
 #     386 "short_than_100nt";
 #     411 "snoRNA";
 #     308 "snRNA";

bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_lv3_shortRNA.bed -wo | wc -l
 # 34
 # two novel miRNAs in the bunch are in gencode

 # chr19   47727524        47727545        novelMiR_556    37.91   -       chr19   47724082        47731482        ENSG00000265134.1       DHS_dyadic      .       miRNA   21
 # chr6    28658285        28658308        novelMiR_210,novelMiR_308,novelMiR_372,novelMiR_593     857.935 -       chr6    28614527        28677863        ENSG00000272278.1       DHS_promoter    .       miRNA   23

bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_lv3_shortRNA.bed -wo > All_miRNA_inter_FANTOM.txt

~/Programs/LiftOverProg/liftOver FANTOM_lv3_shortRNA.bed ~/Programs/LiftOverProg/hg19ToHg38.over.chain "FANTOM_lv3_shortRNA_hg38.bed" "FANTOM_lv3_shortRNA_hg38unmapped.bed"
 # 12 didnt' map 

bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_lv3_shortRNA_hg38.bed -wo > All_miRNA_inter_FANTOM_hg38.txt
```

I am unsure if the "robust" miRNAs in the above FANTOM dataset are included, the same, or different from the atlas of miRNA publication in 2017. I downloaded the supplementary tables with the predicted miRNAs.

make a bed from `FANTOM_miRNA_coordinates.csv`

```R
FANT_coords <- read.csv("FANTOM_miRNA_coordinates.csv")
FANT_coords <- FANT_coords[,c(3,4,5,2,1,6)]
write.table(FANT_coords, "FANTOM_miRNA_coordinates.bed", append=F, col.names=F, row.names=F, quote=F, sep="\t")

FANT_cat <- read.csv("FANTOM_miRNA_set_list.csv")
miRNA_annot <- merge(x=FANT_coords, y=FANT_cat, by="Name")
dim(miRNA_annot)
[1] 6543   14
 # well.. there should have been way more overlap.
table(miRNA_annot$Set)
 candidate permissive     robust 
      6543          0          0 
WHat's the naming scheme for robust?
sum(FANT_cat$Set=="robust") # 795
sum(FANT_cat$Set=="permissive") # 1076
 # ok, that's because the robust and permissive are miRBase annotations

miRNA_annot_permissive <- miRNA_annot[miRNA_annot$Sufficient.tags=="yes",]
nrow(miRNA_annot_permissive)
[1] 321

library(genefilter)
idx_robust <- which(rowSums(miRNA_annot_permissive=="yes") >=4)
length(idx_robust)
 # [1] 282

miRNA_annot_robust <- miRNA_annot_permissive[idx_robust,]
miRNA_annot_permissive <- miRNA_annot_permissive[-idx_robust,]

colnames(miRNA_annot_robust)[8:14] <- c("Suff_tags", "end_consistent", "end_overhang", "suff_bp_nucs", "sufficient_bp_energy", "CAGE_pval_FANTOM", "CAGE_pval_ENCODE")

colnames(miRNA_annot_permissive)[8:14] <- c("Suff_tags", "end_consistent", "end_overhang", "suff_bp_nucs", "sufficient_bp_energy", "CAGE_pval_FANTOM", "CAGE_pval_ENCODE")

write.table(miRNA_annot_robust, "FANTOM_robust_annot.txt", append=F, col.names=F, row.names=F, quote=F, sep="\t")

write.table(miRNA_annot_permissive, "FANTOM_permissive_annot.txt", append=F, col.names=F, row.names=F, quote=F, sep="\t")

miRNA_robust_bed <-miRNA_annot_robust[,c(2,3,4,1,5,6)]
miRNA_permissive_bed <-miRNA_annot_permissive[,c(2,3,4,1,5,6)]

write.table(miRNA_robust_bed, "FANTOM_robust_miRNA_2017.bed", append=F, col.names=F, row.names=F, quote=F, sep="\t")

write.table(miRNA_permissive_bed, "FANTOM_permissive_miRNA_2017.bed", append=F, col.names=F, row.names=F, quote=F, sep="\t")
```

Unfortunately, FANTOM used hg19, so I need to liftover

```bash
~/Programs/LiftOverProg/liftOver FANTOM_robust_miRNA_2017.bed ~/Programs/LiftOverProg/hg19ToHg38.over.chain "FANTOM_robust_miRNA_2017_hg38.bed" "FANTOM_robust_miRNA_2017_hg38unmapped.bed"
 # all lifted
 
~/Programs/LiftOverProg/liftOver FANTOM_permissive_miRNA_2017.bed ~/Programs/LiftOverProg/hg19ToHg38.over.chain "FANTOM_permissive_miRNA_2017_hg38.bed" "FANTOM_permissive_miRNA_2017_hg38unmapped.bed"
 # all lifted
 
 # how many of the 2017 pub are found in the catalog I used before?
bedtools intersect -a FANTOM_permissive_miRNA_2017.bed -b FANTOM_lv3_shortRNA.bed -wo | wc -l # 0

bedtools intersect -a FANTOM_robust_miRNA_2017.bed -b FANTOM_lv3_shortRNA.bed -wo | less wc -l #8

 # not much overlap between nc paper and the miRNA paper
```

Now I can overlap with the miRNAs

```
 # I'm intersecting without enforcing strandedness.. but enforcing it does reduce the amount of overlap.
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_permissive_miRNA_2017_hg38.bed -wo | wc -l
# 3

bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_robust_miRNA_2017_hg38.bed -wo | less

bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_robust_miRNA_2017_hg38.bed -wo > All_miRNA_inter_FANTOM_robust2017.txt

bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/FANTOM_permissive_miRNA_2017_hg38.bed -wo > All_miRNA_inter_FANTOM_permissive2017.txt
```


7. Annotate the miRNAs

downloaded the 5UTR, exon, intron, and 3UTR from the UCSC table browser 
human>hg38>genes and gene predictions > ALl Gencode v24 > Comprehensive 
then hit get out put and selected 5UTR, exon, intron, or 3UTR
based off: https://www.biostars.org/p/94823/

I made 4 bed files in: `~/Genome_annotations_hg38/Gencode_beds`

```bash
# which are intergenic?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/cat_UTR_intron_exon_sort_merge.bed -v | wc -l
# 64, counts chrEBV and JTFH (decoy) miRNAs 
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/cat_UTR_intron_exon_sort_merge.bed -v > All_miRNA_inter_intergenic.txt

 # which are in exons?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_exon_hg38.bed -wo | wc -l
 # 3989
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_exon_hg38.bed -u | wc -l
 # 2856
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_exon_hg38.bed -u > All_miRNA_inter_exon.txt

 # which are in introns?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_intron_hg38.bed -wo | wc -l
 # 10590, so obviously multiple intersects 
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_intron_hg38.bed -u | wc -l
 # 2106
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_intron_hg38.bed -u > All_miRNA_inter_intron.txt

 # which are in 5UTR?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_5UTR_hg38.bed -wo | wc -l
 # 1954, so obviously multiple intersects 
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_5UTR_hg38.bed -u | wc -l
 # 1458
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_5UTR_hg38.bed -u > All_miRNA_inter_UTR5.txt

 # which are in 3UTR?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_3UTR_hg38.bed -wo | wc -l
 # 1847, so obviously multiple intersects 
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_3UTR_hg38.bed -u | wc -l
 # 1536
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_3UTR_hg38.bed -u > All_miRNA_inter_UTR3.txt

 # in the direction of a gene?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_genes_hg38.bed -s -u | wc -l 
 # 2877
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_genes_hg38.bed -s -u > All_miRNA_inter_gene_stranded.txt

 # how many overlap gene annotation?
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_genes_hg38.bed -u | wc -l 
 # 2970
bedtools intersect -a All_miRNA_chrfix.bed -b ~/Genome_annotations_hg38/Gencode_beds/v24_all_genes_hg38.bed -u > All_miRNA_inter_gene.txt

 # Annotate the miRNA by ATACseq peaks
bedtools intersect -a All_miRNA_chrfix.bed -b ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_numPeak.bed -wo | wc -l 
 # 440, chrEBV miRNAs are included 
bedtools intersect -a All_miRNA_chrfix.bed -b ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_numPeak.bed -wo > All_miRNA_inter_ATACpeak.txt

 # How many are located +/-250bp from a summit
bedtools intersect -a All_miRNA_chrfix.bed -b ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_summit_noNeg.bed -wo | wc -l 
 # 199 miRNAs are in summits
bedtools intersect -a All_miRNA_chrfix.bed -b ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_summit_noNeg.bed -wo > All_miRNA_inter_ATACsummit.txt

 # does a peak precedes the miRNA and in the right direction.
module load bedtools2
module load bedops/2.4.12/gcc.4.9.2

# Need to sort the ATACpeaks and the miRNAs
bedtools sort -i ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_numPeak.bed > ATACpeaks_sorted.bed
bedtools sort -i All_miRNA_chrfix.bed > All_miRNA_chrfix_sorted.bed

 # find the closest feature
closest-features All_miRNA_chrfix_sorted.bed ATACpeaks_sorted.bed | less 

closest-features All_miRNA_chrfix_sorted.bed ATACpeaks_sorted.bed | awk '{if ($6~"-" && ($12-$3) < 2000) print $4"\t"$14"\t"$12-$3; if ($6~"+" && ($2-$8) < 2000) print $4"\t"$9"\t"$2-$8}' | awk 'NR>2' > All_miRNAs_with_ATAC_annot.txt 

wc -l All_miRNAs_with_ATAC_annot.txt
 # 619
```

 # getting the mature sequence for the final set of novel miRNAs 
for pooled in either: (All_novel_chr16.bed, All_novel_chr17_1half.bed, All_novel_chr17_2half.bed, All_novel_chr18up.bed, All_novel_chr21_1half.bed, All_novel_chr21_2half.bed)
all.filtered.primary.chr17.1half.result 
all.filtered.primary.chr17.2half.minusMIR21.result
all.filtered.primary.chr18up.result 

```bash
 # novelMiR_373
grep "novelMiR_373" all.filtered.primary.chr18up.result | less 

 # novelMiR_268 - found from both pooled and single cell (I used the single individual for naming.. but unsure which cell line)
grep "novelMiR_268"  GM93Aligned.sortedByCoord.out.filtered.result | less

 # novelMiR_19 
grep "novelMiR_19" all.filtered.primary.chr17.2half.minusMIR21.result  | less 

 #novelMiR_383
grep "novelMiR_383"  GM80Aligned.sortedByCoord.out.filtered.result | less

 # novelMiR_550
grep "novelMiR_550" all.filtered.primary.chr16.result| less 

 # novelMiR_948
grep "novelMiR_948" all.filtered.primary.chr16.result| less 

 # novelMiR_523
grep "novelMiR_523" all.filtered.primary.chr16.result| less 

 # novelMiR_86_pooled
grep "novelMiR_86" all.filtered.primary.chr18up.result | less 

 # novelMiR_2367
grep "novelMiR_2367"  all.filtered.primary.chr16.result| less 

 # novelMiR_2083
grep "novelMiR_2083"  all.filtered.primary.chr16.result| less

 # novelMiR_862
grep "novelMiR_862"  all.filtered.primary.chr16.result| less  

 # novelMiR_2667
grep "novelMiR_2667"  all.filtered.primary.chr16.result| less 

 # novelMiR_2367
grep "novelMiR_2367" all.filtered.primary.chr16.result| less 

 ### #####
 # ENSG00000206625.1 = U6

grep "ENSG00000206625.1" GM77_refSeq_filt.txt 
 # roughly correlates with the depth of coverage (unsure if it correlates with miRNA coverage)

less /home/greally-lab/indexes/hg38/miRNA/mature.hsa.dna.fa

grep -e "MIMAT0002174" -e "MIMAT0000693" -e "MIMAT0004518" -e "MIMAT0004672" -e "MIMAT0005951" /home/greally-lab/indexes/hg38/miRNA/hsa.gff3 

grep -e "hsa-miR-1307-3p" -e "hsa-miR-106b-3p" -e "hsa-miR-16-2-3p" -e "hsa-miR-30e-3p" -e "hsa-miR-484" -A 1 /home/greally-lab/indexes/hg38/miRNA/mature.hsa.dna.fa

grep "miR-155" /home/greally-lab/indexes/hg38/miRNA/hsa.gff3 
```