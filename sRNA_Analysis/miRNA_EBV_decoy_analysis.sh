## Redoing miRNA analysis with the EBV_decoy GENCODE!!! genome 
cat GRCh38.primary_assembly.genome.fa.gz GCA_000786075.2_hs38d1_genomic.fna.gz chrEBV.fa.gz > Gencode_EBV_decoy.fa.gz
mv Gencode_EBV_decoy.fa.gz ../hg38_GenCode_EBV_decoy/
gunzip -d Gencode_EBV_decoy.fa.gz
# pwd = /home/greally-lab/indexes/hg38_GenCode_EBV_decoy
## Generating the EBV_decoy genome for miRDeep
sed -e 's/^\(>[^[:space:]]*\).*/\1/' Gencode_EBV_decoy.fa > Gencode_EBV_decoy_noComments.fa
mkdir hg38_GenCode_EBV_decoy
mv Gencode_EBV_decoy_noComments.fa hg38_GenCode_EBV_decoy/
perl /home/greally-lab/indexes/Toxo_Human/split.pl Gencode_EBV_decoy_noComments.fa

# pwd /home/ajohnsto/Programs/MDS_command_line_v37/MDS_command_line/
mv /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/hg38_GenCode_EBV_decoy genome/

qsub -S /bin/bash -N decoy_build -cwd -q highmem.q -l h_vmem=100G -j y -pe smp-highmem 1 << EOF
module load MDS
java -jar -Xmx75G build_bwt_idx.jar hg38_GenCode_EBV_decoy
EOF


## Align to the EBV_decoy genome
for f1 in *_comb.fq.gz
do
FILE="$(echo ${f1} | cut -d '/' -f2)"
SAMPLE="$(echo ${f1} | cut -d '_' -f1)"	
echo ${SAMPLE}
qsub -S /bin/bash -N ${SAMPLE}_alignGen_decoy -cwd -q highmem.q -l h_vmem=5.6G -j y -pe smp-highmem 23 << EOF
module load samtools
module load STAR
STAR --runThreadN 23 --genomeDir /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Star --readFilesIn ${f1} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 --outFileNamePrefix Mapped_GenCode_EBV_decoy/${SAMPLE}
EOF
done

## taking out reads that were soft clipped in the 5' end
qrsh -q highmem.q -l h_vmem=20G -pe smp-highmem 1
for f1 in *Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
samtools view -h ${f1} | awk '{S=0; split($6,C,/[0-9]*/); n=split($6,L,/[NMSID]/);  if (and($2,0x10)>0 && C[n]=="S") {S=L[n-1]} else if (and($2,0x10)==0 && C[2]=="S") {S=L[1]}; if (S<=1) print }' | samtools view -Sb - > ${SAMPLE}Aligned.sortedByCoord.out.filtered.bam
done

## Run the filtered reads through miRDeep*
for f1 in ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy/*filtered.bam
do
dest1="${f1##*/}"
SAMPLE="$(echo ${dest1} | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N ${SAMPLE}_miRDeep -cwd -q highmem.q -l h_vmem=100G -j y -pe smp-highmem 1 << EOF
module load MDS
module load samtools
java -jar -Xmx75g MD.jar -g hg38_GenCode_EBV_decoy -t 16 -l 25 ${f1}
EOF
done

## Count refseq GTF using primary alignment
for f1 in *Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N ${SAMPLE}_featureCounts -j y -cwd -q highmem.q -pe smp-highmem 4 -l h_vmem=10G << EOF
module load subread/1.5.0-p1/gcc.4.4.7 
featureCounts  -T 4 -F GTF --primary -a /home/greally-lab/indexes/hg38_EBV/gencode.v24.primary_assembly.annotation.EBV.gtf -o Count_refseq_primary/${SAMPLE}_refSeq.txt ${f1}
EOF
done

## Count rRNA/tRNA using primary alignment
for f1 in *Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N ${SAMPLE}_featureCounts -j y -cwd -q highmem.q -pe smp-highmem 4 -l h_vmem=10G << EOF
module load subread/1.5.0-p1/gcc.4.4.7 
featureCounts -T 4 -F GTF --primary -a /home/greally-lab/indexes/hg38/rRNA/rRNAandtRNA.GTF -o Count_rRNA_n_tRNA_primary/${SAMPLE}_rRNA.txt ${f1}
EOF
done

## Discover novel miRNAs using miRDeep
# Find number of known miRNAs
touch Known_found.txt 
for f1 in *Aligned.sortedByCoord.out.filtered.result;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
awk -v x=$SAMPLE '$1~"hsa"{count++}END{print x"\t"count}' ${f1} >> Known_found.txt 
done 

# find number of novel
touch Novel_found.txt 
for f1 in *Aligned.sortedByCoord.out.filtered.result;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
awk -v x=$SAMPLE '$1~"novel"{count++}END{print x"\t"count}' ${f1} >> Novel_found.txt 
done 

## determine a cutoff for novel transcripts
qrsh -l h_vmem=5G
# what are the means of novel vs known? 
grep hsa GM88Aligned.sortedByCoord.out.filtered.result | awk '{ sum += $2 } END { if (NR > 0) print sum / NR }'
#11540.6
grep novel GM88Aligned.sortedByCoord.out.filtered.result | awk '{ sum += $2 } END { if (NR > 0) print sum / NR }'
#1413.56

# the mean is heavily skewed by a couple miRNAs so, will look at median

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

# finding median
awk '$1~"hsa"{print $2}' GM88Aligned.sortedByCoord.out.filtered.result | sort -n | awk -f median.awk 
# 2.11
awk '$1~"novel" {print $2}' GM88Aligned.sortedByCoord.out.filtered.result | sort -n | awk -f median.awk 
#-2.85

# number of novel
grep novel GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
#976
# number of known
grep hsa GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
#467

# Odds score > x for known
awk '$1~"hsa" && $2>10 {print $2}' GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
#201
awk '$1~"hsa" && $2>0 {print $2}' GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
# 341

# Odds score > x for novel
awk '$1~"novel" && $2>10 {print $2}' GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
# 61
awk '$1~"novel" && $2>0 {print $2}' GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
# 238

awk '$1~"novel" && $2>0 {print $6}' GM88Aligned.sortedByCoord.out.filtered.result | sort -n | less

awk '$1~"novel" && $6>50 && $2>0 {print $6}' GM88Aligned.sortedByCoord.out.filtered.result | wc -l 
# 223

touch Novel_filtered.txt 
for f1 in *Aligned.sortedByCoord.out.filtered.result;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
awk -v x=$SAMPLE '$1~"novel" && $6>50 && $2>0{count++}END{print x"\t"count}' ${f1} >> Novel_filtered.txt 
done 

## merge the novel transcripts from all 17 samples
module load bedtools2 

# Need to create a bed file out of the filtered transcripts then need to look at intersections
awk '$1~"novel" && $6>50 && $2>0 {print $7}' GM88Aligned.sortedByCoord.out.filtered.result | awk -F "-" '{print $1"\t"$2}' > TestStart_stop.txt
awk '$1~"novel" && $6>50 && $2>0 {print $3"\t"$1"\t"$2"\t"$4}' GM88Aligned.sortedByCoord.out.filtered.result | paste - TestStart_stop.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > Test.bed

for f1 in *Aligned.sortedByCoord.out.filtered.result;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
awk '$1~"novel" && $6>50 && $2>0 {print $7}' ${f1} | awk -F "-" '{print $1"\t"$2}' > ${SAMPLE}_start_stop.txt
awk '$1~"novel" && $6>50 && $2>0 {print $3"\t"$1"\t"$2"\t"$4}' ${f1} | paste - ${SAMPLE}_start_stop.txt | awk -F "\t" '{print $1"\t"$5"\t"$6"\t"$2"\t"$3"\t"$4}' > ${SAMPLE}_novel_filtered.bed 
done

# intersecting the bed files 
bedtools intersect -wo -a GM83_novel_filtered.bed -b GM80_novel_filtered.bed | less
# 13th column has the intersection bp overlap. most are 20 bp and up
cat *_novel_filtered.bed > Novel_filt_bed_list.bed
sortBed -i Novel_filt_bed_list.bed >  Novel_filt_bed_list_sort.bed

bedtools merge -s -i Novel_filt_bed_list_sort.bed > Merged_novel_fillist.bed
#139
bedtools merge -s -c 6,5,4 -o distinct,mean,distinct -i Novel_filt_bed_list_sort.bed > Merged_novel_fillist.bed
#139

## Make new GTF file containing novel transcripts
# 1 repeat that I needed to manually fix
awk '{print $6}' Merged_novel_fillist.bed | awk -F "," '{print $1}' | sort | uniq -d 
vi Merged_novel_fillist.bed
awk '{print $6}' Merged_novel_fillist.bed | awk -F "," '{print $1}' |  awk '{sub($1, "\"&\""); print}' > Merged_novel_fillist_names.txt


awk '{print $1"\t.\tmiRNA\t"$2"\t"$3"\t"$4"\t"$5"\t.\tID"}' Merged_novel_fillist.bed > Merged_novel_fillist_gtf.txt

paste -d " " Merged_novel_fillist_gtf.txt Merged_novel_fillist_names.txt  > Merged_novel_fillist.gtf

## need to remove the chr off non chromosome contigs
# had to manually take out the chr in front of the non-chr contigs such as GL000220.1 
cp Merged_novel_fillist.gtf Merged_novel_fillist_nochr.gtf

vi Merged_novel_fillist_nochr.gtf

cat /home/greally-lab/indexes/hg38/miRNA/hsa.GTF Merged_novel_fillist_nochr.gtf > miRNA_hg38_EBV_decoy.gtf

sort miRNA_hg38_EBV_decoy.gtf > miRNA_hg38_EBV_decoy_sort.gtf

wc -l miRNA_hg38_EBV_decoy_sort.gtf
# 4833 miRNA_hg38_EBV_decoy_sort.gtf
# only 2952 that are miRNAs

## count the miRNAs using featureCounts - not using primary because I want the most stringent counting for the miRNAs. I will use the unfiltered bam file because I'm more curious about counts not the eqxact positioning of the miRNA transcript

for f1 in *Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N ${SAMPLE}_featureCounts -j y -cwd -pe smp 2 -l h_vmem=5.6G << EOF
module load subread/1.5.0-p1/gcc.4.4.7 
featureCounts -T 2 -F GTF -t miRNA -g ID -a miRNA_hg38_EBV_decoy_sort.gtf -o Count_miRNA/${SAMPLE}_miRNA.txt ${f1}
EOF
done

awk 'NR>2 {print $2}' GM77_miRNA.txt | uniq -c

awk '{print $3}' miRNA_hg38_EBV_decoy_sort.gtf | sort | uniq -c
# 2952 miRNA
# 1881 miRNA_primary_transcript

 qstat -u "*" | awk '{print $4}' | sort | uniq -c

17,808,581
8,000,000
61,154,446
