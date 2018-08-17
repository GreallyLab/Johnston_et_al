# testing if pre-filtering affects the results of miRDeep* because I had a large number of reads with secondary alignments, which could have led to many "novel" miRNAs being the result of mis-mapping

# miRdeep sorts and indexes the bam files first. Let's flagstat to find out. gm79 has least amount of reads and is the best test sample

qsub -S /bin/bash -N Test_flagstat -j y -cwd -l h_vmem=10G << EOF
module load samtools
samtools flagstat GM85Aligned.sortedByCoord.out.filtered.bam > GM85Aligned.sortedByCoord.out.filtered.bam.flagstat
samtools flagstat GM85Aligned.sortedByCoord.out.filtered.sorted.bam > GM85Aligned.sortedByCoord.out.filtered.sorted.bam.flagstat
EOF 

# GM85Aligned.sortedByCoord.out.filtered.bam
183,378,711 + 4,798,420 in total (QC-passed reads + QC-failed reads)
144914425 + 3744914 secondary
0 + 0 supplementary
0 + 0 duplicates
183378711 + 4798420 mapped (100.00%:100.00%)

# GM85Aligned.sortedByCoord.out.filtered.sorted.bam.flagstat
48,709,774 + 1485506 in total (QC-passed reads + QC-failed reads)
33,186,282 + 1035316 secondary
0 + 0 supplementary
0 + 0 duplicates
48709774 + 1485506 mapped (100.00%:100.00%)

# So miRdeep doesn't remove secondary alignments. It may take all reads that are best hit though since the number of reads is now just a little more than 48,006,

# keep only the uniquely aligned reads 

for f1 in GM*Aligned.sortedByCoord.out.filtered.bam
do
SampleName="$(echo ${f1} | cut -c 3-4)"
echo $SampleName
qsub -S /bin/bash -N filter_${SampleName}_StarBam -l h_vmem=10G -j y -cwd << EOF
module load samtools/1.2/gcc.4.4.7
# Filtered unmapped reads, not primary alignment, not passing vendor check, and duplicates
samtools view -b -F 512 -F 4 -F 256 ${f1} > ${f1%.bam}.primary.bam
EOF
done

# Then merge the uniquely aligned reads 
qsub -S /bin/bash -N Combine_primaryBam -j y -cwd -l h_vmem=50G << EOF
module load samtools
samtools merge All.Aligned.sortedByCoord.out.filtered.primary.bam *Aligned.sortedByCoord.out.filtered.primary.bam
EOF
# 17G Jan 27 06:48 All.Aligned.sortedByCoord.out.filtered.primary.bam
qsub -S /bin/bash -N flagstat_primaryBam -j y -cwd -l h_vmem=10G << EOF
module load samtools
samtools flagstat All.Aligned.sortedByCoord.out.filtered.primary.bam >  All.Aligned.sortedByCoord.out.filtered.primary.flagstat
EOF

# 657,893,266 million reads are primary mapped. 

qsub -S /bin/bash -N all_primary_miRDeep -q all.q@@r2940.v4 -cwd -l h_vmem=240G -l exclusive=true -j y << EOF
module load MDS
module load samtools
java -jar -Xmx210g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/All.Aligned.sortedByCoord.out.filtered.primary.bam
EOF

# I ran out of memory using this
# all.q@n902 is the new highmem that I should try

qrsh -q all.q@@n902  -l h_vmem=1000G

qrsh -q highmem.q -l h_vmem=1430G


qsub -S /bin/bash -N all_primary_miRDeep -q highmem.q -l h_vmem=1430G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/All.Aligned.sortedByCoord.out.filtered.primary.bam
EOF

# this is still not enough memory :( 
# So let's try to split up the chromosomes again!

module load samtools 
samtools view All.Aligned.sortedByCoord.out.filtered.primary.bam 1:1-999999999 | less

module load bamtools
bamtools split -h
bamtools split -in All.Aligned.sortedByCoord.out.filtered.primary.bam -reference


samtools merge all.filtered.primary.chr16.bam *.bam

qrsh -l h_vmem=15G

samtools merge all.filtered.primary.chr17up.bam *.bam
# need to make a 17up canon only + EBV +chr2
mv *chr* ../split_primary_can17/

# I needed to further split the bam file

# I took out chr21 as it was the largest then combined 17up+2+EBV - 21
qsub -S /bin/bash -N merge_bam -q highmem.q -l h_vmem=50G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.chr17up.bam *.bam
EOF 

# this was too much for the computer. ALSO, it looks like I get the same results for each miRNA regardless if chromosome is analyzed seperately or together. 

# took out chr 17, 21.. now its 18 and up
qsub -S /bin/bash -N merge_bam -l h_vmem=6G -cwd -j y << EOF
module load samtools
samtools merge all.filtered.primary.chr18up.bam *.bam
EOF 

# chr17 and 21 seem to not be working by themselves, so I'm going to split them in 2.
# chr 21 has a length of  46709983 , so half is 23354990
qsub -S /bin/bash -N index_bam -l h_vmem=6G -cwd -j y << EOF
module load samtools
samtools index All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.bam
EOF
samtools view -b All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.bam chr21:1-23354990 > All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.bam
# chr21 -  second half
samtools view -b All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.bam chr21:23354990-46709983 > All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.bam
# chr21 -  second half is still causing memory issues.. split 1 more time.
# 36,017,442 reads
samtools index All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.bam
awk '$1=="chr21" && $4>23354990 && $5<46709983' hsa.GTF | less

# MIR 155 
samtools view All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.bam chr21:25573980-25574044 | wc -l # 35,589,933 reads (which is the issue)

samtools view -b All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.bam chr21:23354990-25573980 chr21:25574044-46709983 > All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.bam


# chr 17 has a length of  83257441 , so half is 41628720
qsub -S /bin/bash -N index_bam -l h_vmem=6G -cwd -j y << EOF
module load samtools
samtools index all.filtered.primary.chr17.bam
EOF
samtools view -b all.filtered.primary.chr17.bam chr17:1-41628720 > all.filtered.primary.chr17.1half.bam
# the second half is still causing memory issues.. split 1 more time.
qsub -S /bin/bash -N index_bam -l h_vmem=3G -cwd -j y << EOF
module load samtools
samtools index all.filtered.primary.chr17.2half.bam
EOF
# 62443080
qsub -S /bin/bash -N split_bam -l h_vmem=2G -cwd -j y << EOF
module load samtools
samtools view -b all.filtered.primary.chr17.2half.bam chr17:41628720-59841266 chr17:59841337-83257441 > all.filtered.primary.chr17.2half.minusMIR21.bam
EOF
samtools view -b all.filtered.primary.chr17.2half.bam chr17:50000000-83257441 > all.filtered.primary.chr17.22half.bam

chr17:59,841,266-59,841,337

awk '$1=="chr17" && $4>59500000 && $5<62443080' hsa.GTF | less
samtools view all.filtered.primary.chr17.2half.bam chr17:59841266-59841337 | wc -l
# 32,373,830 million which is microRNA 21.. which is breaking the program. 
samtools flagstat all.filtered.primary.chr17.2half.bam # 37,653,119

samtools index all.filtered.primary.chr17.2half.minusMIR21.bam
all.filtered.primary.chr17.2half.minusMIR21.bam
samtools view all.filtered.primary.chr17.2half.minusMIR21.bam chr17:59841266-59841337 | wc -l

samtools view all.filtered.primary.chr17.2half.minusMIR21.bam | wc -l
# 5,285,593


# from chr 16 down
qsub -S /bin/bash -N all_primary_chr16down -q highmem.q -l h_vmem=1430G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can/all.filtered.primary.chr16.bam
EOF

#chr17 - 1st half
qsub -S /bin/bash -N all_primary_chr171half -q highmem.q -l h_vmem=1430G -l exclusive=true -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can17/all.filtered.primary.chr17.1half.bam
EOF

# chr17 - 2nd half
qsub -S /bin/bash -N all_primary_chr17_2half -q highmem.q -l h_vmem=1430G -l exclusive=true -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can17/all.filtered.primary.chr17.2half.bam
EOF

# chr17 - 2nd half minus MIR 21
qsub -S /bin/bash -N all_primary_chr17_2half_minus -q highmem.q -l h_vmem=1430G -l exclusive=true -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can17/all.filtered.primary.chr17.2half.minusMIR21.bam
EOF

# chr21 - 1st half
qsub -S /bin/bash -N all_primary_chr21_1half -q highmem.q -l exclusive=true -l h_vmem=1430G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can21/All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.1half.bam
EOF

# chr21 - 2nd half minus MIR 155
qsub -S /bin/bash -N all_primary_chr21_2half_minus -q highmem.q -l exclusive=true -l h_vmem=1430G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can21/All.Aligned.sortedByCoord.out.filtered.primary.REF_chr21.2half.minusMIR155.bam
EOF

# 18up minus 21 # This worked
qsub -S /bin/bash -N all_primary_chr21 -q highmem.q -l exclusive=true -l h_vmem=1430G -cwd -j y << EOF
module load MDS
module load samtools
java -jar -Xmx1350g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/split_primary_can18up/all.filtered.primary.chr18up.bam
EOF


## median.awk
#/usr/bin/env awk
{
    count[NR] = $1;
}
END {
    if (NR % 2) {
        print count[(NR + 1) / 2];
    } else {
        print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;
    }
}
###


wc -l All.Aligned.sortedByCoord.out.filtered.primary.result
# 1917 All.Aligned.sortedByCoord.out.filtered.primary.result

grep novel All.Aligned.sortedByCoord.out.filtered.primary.result| wc -l 
# 1653

# number of known
grep hsa All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l 
# 264

awk '$1~"hsa"{print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | sort -n | awk -f median.awk 
# 3.785 (15.33)
awk '$1~"novel" {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | sort -n | awk -f median.awk 
# -1.60 (-2.08)


awk '$1~"hsa" && $2>10 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l 
# 127 (151)
awk '$1~"hsa" && $2>0 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l 
# 195 (219)

awk '$1~"novel" && $2>10 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l 
# 191 (246)
awk '$1~"novel" && $2>0 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l 
# 568 (704)

awk '$1~"novel" && $6>50 && $2>10 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l #59 novel miRNA

awk '$1~"novel" && $6>50 && $2>0 {print $2}' All.Aligned.sortedByCoord.out.filtered.primary.result | wc -l #126 novel miRNA

awk '$1~"novel" && $6>50 && $2>10 {print $7}' All.Aligned.sortedByCoord.out.filtered.primary.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop10.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' All.Aligned.sortedByCoord.out.filtered.primary.result | paste - AllStart_stop10.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel10.bed # 59



awk '{print $3}' All.Aligned.sortedByCoord.out.filtered.primary.result | sort | uniq -c

awk '{print $3}' all.filtered.primary.chr17up.result | sort | uniq -c
# 107 chr17
# 172 chr2

awk '{print $3}' all.filtered.primary.chr16.result | sort | uniq -c

wc -l All_novel10.bed
# 59

1 chr
    228 chr1
    103 chr10
    128 chr11
    115 chr12
     31 chr13
    106 chr14
     72 chr15
    136 chr16
    107 chr17
    172 chr2
    131 chr3
     68 chr4
    115 chr5
    119 chr6
    106 chr7
     82 chr8
     98 chr9

awk '$1~"novel" && $6>50 && $2>0 {print $7}' All.Aligned.sortedByCoord.out.filtered.primary.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop0.txt

awk '$1~"novel" && $6>50 && $2>0 {print $3"\t"$1"\t"$2"\t"$4}' All.Aligned.sortedByCoord.out.filtered.primary.result | paste - AllStart_stop0.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_primary_0.bed # 126

module load bedtools2
# intersect with the novel miRNAs found before
bedtools intersect -wo -a All_novel.bed -b Merged_novel_fillist.bed | wc -l
# only 23 :-/.. not great overlap

bedtools intersect -wo -a All_novel.bed -b novel_miRNA_final.bed | wc -l
# 3

bedtools intersect -wo -a novel_miRNA_final.bed -b All_novel10.bed | less
# only one overlapping :( 
# chr14   49586857        49586879        novelMiR_103,novelMiR_124,novelMiR_174,novelMiR_203,novelMiR_247,novelMiR_68,novelMiR_95        1383.38 +       chr14   49586857        49586879        novelMiR_1322   5824.83 +       22
bedtools intersect -wo -a novel_miRNA_final.bed -b All_novel_primary_0.bed | less
# only 2.. 
# chr12   80832562        80832583        novelMiR_157,novelMiR_193,novelMiR_235  17.95333333     +       chr12   80832562        80832583        novelMiR_1238   1.90    +       21

 
awk '{print $1}' novel_miRNA_final.bed | sort | uniq -c

awk '{print $3}' all.filtered.primary.chr17up.result | sort | uniq -c

all.filtered.primary.chr17up.result

wc -l all.filtered.primary.chr16.result
# 3212 

grep novel all.filtered.primary.chr16.result | wc -l 
# 2748

# number of known
grep hsa all.filtered.primary.chr16.result | wc -l 
# 463

awk '$1~"hsa"{print $2}' all.filtered.primary.chr16.result | sort -n | awk -f median.awk 
# 2.12 (3.785) (15.33)
awk '$1~"novel" {print $2}' all.filtered.primary.chr16.result | sort -n | awk -f median.awk 
# -1.675 (-1.60) (-2.08)

awk '$1~"novel" && $6>50 && $2>10 {print $2}' all.filtered.primary.chr16.result | wc -l #104 (59) novel miRNA

awk '$1~"novel" && $6>50 && $2>0 {print $2}' all.filtered.primary.chr16.result| wc -l #220 (126) novel miRNA


awk '$1~"novel" && $6>50 && $2>10 {print $7}' all.filtered.primary.chr16.result | awk -F "-" '{print $1"\t"$2}' > AllStart_stop_chr16.txt

awk '$1~"novel" && $6>50 && $2>10 {print $3"\t"$1"\t"$2"\t"$4}' all.filtered.primary.chr16.result | paste - AllStart_stop_chr16.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > All_novel_chr16.bed # 59


bedtools intersect -wo -a ../novel_miRNA_final.bed -b All_novel_chr16.bed | less
bedtools intersect -wo -a ../novel_miRNA_final.bed -b All_novel_chr16.bed | wc -l
#10


novel_miRNA_final.bed

awk '{print ($3-$2)}' novel_miRNA_final.bed


awk '{count += ($3-$2)} END {print count/NR}' novel_miRNA_final.bed # 20.9167 mean length




