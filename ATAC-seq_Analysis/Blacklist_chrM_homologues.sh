# using wgsim to simulate mitochondrial reads to then map to the genome and create a peak blacklist (this blacklist is created because there are mitochondrial homologues in the reference genome)

# getting the mitochondrial chromosome 
# in  /home/greally-lab/indexes/hg38_GenCode_EBV_decoy

module load samtools 
samtools faidx Gencode_EBV_decoy_noComments.fasta chrM > chrM.fa

cat Gencode_EBV_decoy_noComments.fasta | awk '{if (substr($0,1) == ">chrM") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' > Gencode_EBV_decoy_noComments_nochrM.fasta

wc -l Gencode_EBV_decoy_noComments_nochrM.fasta
51,740,953 Gencode_EBV_decoy_noComments_nochrM.fasta
wc -l chrM.fa
278 chrM.fa
wc -l Gencode_EBV_decoy_noComments.fasta
51,741,231 Gencode_EBV_decoy_noComments.fasta

# checks out! 

# need to make a bwa index that doesn't include chrM.
qsub -S /bin/bash -N Build_ref_nochrM -cwd -l h_vmem=30G -j y << EOF
module load bwa 
bwa index -a bwtsw Gencode_EBV_decoy_noComments_nochrM.fasta 
EOF

# in anjohnst/Programs/wgsim-master

# https://github.com/lh3/wgsim/

./wgsim -d 400 -e 0 -N 10000000 -1 100 -2 100 -R 0 -r 0 /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/chrM.fa chrM.read1.fq chrM.read2.fq

PICARD_ROOT='/public/apps/picard-tools/1.92/'
module load picard-tools/1.92/java.1.8.0_20 
qsub -S /bin/bash -N chrM_align -j y -cwd -pe smp 20 -l h_vmem=5.6G << EOF
module load bwa/0.7.13/gcc.4.4.7
bwa mem -M -t 20 /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments_nochrM.fasta  chrM.read1.fq chrM.read2.fq > chrM.sam

module load samtools
module load picard-tools/1.92/java.1.8.0_20 
module load MACS2/2.1.0/python.2.7.8
module load bedtools2
module load samtools 

## Select only uniquely mapped reads and sort
samtools view -bq 1 -F 4 chrM.sam | samtools sort -@ 15 - chrM_UMap
echo "chrM_UMap.bam" >> chrM_flagstat.txt
samtools flagstat chrM_UMap.bam >> chrM_flagstat.txt
samtools index chrM_UMap.bam
	
## Remove duplicates
java -jar $(which MarkDuplicates.jar) INPUT=chrM_UMap.bam OUTPUT=chrM_UMap_mkdup.bam METRICS_FILE=chrM_UMap_noMT_mkdup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates ASSUME_SORTED=true MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000
echo "chrM_UMap_mkdup.bam" >> chrM_flagstat.txt
samtools flagstat chrM_UMap_mkdup.bam >> chrM_flagstat.txt
samtools index chrM_UMap_mkdup.bam

samtools view -h -f 0x0040 -b chrM_UMap_mkdup.bam > chrM_UMap_mkdup_read1.bam
bedtools bamtobed -i chrM_UMap_mkdup_read1.bam > chrM_UMap_mkdup_read1.bed
awk 'BEGIN {OFS = "\t"} ; {if (\$6 == "+") print \$1, \$2 + 4, \$3 + 4, \$4, \$5, \$6; else print \$1, \$2 - 5, \$3 - 5, \$4, \$5, \$6}' chrM_UMap_mkdup_read1.bed > chrM_UMap_mkdup_read1_shifted.bed

#Jason's method
macs2 callpeak --nomodel -t chrM_UMap_mkdup_read1_shifted.bed -n chrM --nolambda -g 3e9 --keep-dup 'all'  
EOF

wc -l chrM_peaks.narrowPeak
# 72 chrM_peaks.narrowPeak vs. 111 peaks found in hg 19 by Greenleaf

# make narrowPeak into a bed file for bedtools intersecting
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' chrM_peaks.narrowPeak > Mito_peak_blacklist.bed





