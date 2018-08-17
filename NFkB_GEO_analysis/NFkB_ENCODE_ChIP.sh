# aligning the NFkB data 



mv ENCFF002EKN.fastq.gz	rep1_2.fq.gz
# 5bd7e13ff8ba9c9551d9bf52c7588cb4 mine
# 5bd7e13ff8ba9c9551d9bf52c7588cb4

mv ENCFF002EKP.fastq.gz rep1_1.fq.gz
# 079969b3024940d10fc2915c87a5e190 mine
# 079969b3024940d10fc2915c87a5e190

mv ENCFF002EKQ.fastq.gz rep2_2.fq.gz
# 29af3314af6a5f385ed255ff839fae39 mine 
# 29af3314af6a5f385ed255ff839fae39

mv ENCFF002EKS.fastq.gz rep2_1.fq.gz 
# ba5a7f3487172b60ba8d5cf0aa24d914 mine
# ba5a7f3487172b60ba8d5cf0aa24d914

ENCFF002EKQ controlled by ENCFF002EFS (paried with ENCFF002EFU) ENCLB793DEG DEG
ENCFF002EKP controlled by ENCFF002EFT
ENCFF002EKS controlled by ENCFF002EFU
ENCFF002EKN controlled by ENCFF002EFQ (paired with ENCFF002EFT)


/Users/Andrew/Documents/Greally\ Lab/sratoolkit.2.8.1-mac64/bin/fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip SRR1171946


module load trim_galore 

qsub -S /bin/bash -N Rep1_trim -j y -cwd -q highmem.q -l h_vmem=20G << EOF
module load trim_galore 
trim_galore --paired --fastqc rep1_1.fq.gz rep1_2.fq.gz
EOF

qsub -S /bin/bash -N Rep2_trim -j y -cwd -q highmem.q -l h_vmem=20G << EOF
module load trim_galore 
trim_galore --paired --fastqc rep2_1.fq.gz rep2_2.fq.gz
EOF

MONK:350:C3TRNACXX:6:2106:3847:75402 1:N:0:CGATGT


# To test phred score using first couple of lines in fastq 
cat test  | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding!";}'

mkdir Bams
mkdir Sams
mkdir Flagstats
mkdir Picard_metrics
module load picard-tools/1.92/java.1.8.0_20 
for f1 in rep*_*val_*.fq.gz
do 
ReadNum="$(echo ${f1} | cut -d '_' -f2)"
echo $ReadNum
SampleName="$(echo ${f1} | cut -d '_' -f1)"
echo $SampleName
if [ $ReadNum != '2' ]
then
Read1="$(echo ${f1})"
echo $Read1	
else
echo ${f1}
qsub -S /bin/bash -N ${SampleName}_align -R y -j y -cwd -q highmem.q -l h_vmem=25G << EOF
module load bwa/0.7.15/gcc.4.4.7
#bwa mem -M -t 1 /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa $Read1 ${f1} > Sams/${SampleName}.sam
module load picard-tools/1.92/java.1.8.0_20
module load samtools/0.1.19/gcc.4.4.7 

# get unique reads 
samtools view -Sbq 1 -F 4 Sams/${SampleName}.sam | samtools sort -@ 1 - Bams/${SampleName}_UMap
echo "${f1%.sam}_UMap.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap.bam >> Flagstats/${SampleName}_flagstat.txt
#samtools index Bams/${SampleName}_UMap.bam

## mark duplicates 
java -jar $(which MarkDuplicates.jar) INPUT=Bams/${SampleName}_UMap.bam OUTPUT=Bams/${SampleName}_UMap_mkdup.bam METRICS_FILE=Picard_metrics/${SampleName}_UMap_noMT_mkdup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates ASSUME_SORTED=true MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000
echo "${SampleName}_UMap_mkdup.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap_mkdup.bam >> Flagstats/${SampleName}_flagstat.txt

# also remove duplicates with samtools rmdup
samtools rmdup Bams/${SampleName}_UMap_mkdup.bam Bams/${SampleName}_UMap_mkdup_rmdup.bam
echo "${SampleName}_UMap_mkdup_rmdup.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap_mkdup_rmdup.bam >> Flagstats/${SampleName}_flagstat.txt
samtools index Bams/${SampleName}_UMap_mkdup_rmdup.bam
EOF
fi
done

# https://github.com/lh3/bwa/issues/102 problem with multithreading

mkdir Peaks
mkdir Summits
cd Bams
for s in *_UMap_noMT_mkdup_rmdup.bam;
#for s in gm83_UMap_noMT_mkdup.bam;
do
SampleName="$(echo ${s} | cut -d '_' -f1)"
DIR='../Peaks/'
DIR2='../Summits/'
out_peaks="$DIR${SampleName}"
out_summits="$DIR2${SampleName}"
qsub -S /bin/bash -N Rep2_trim -j y -cwd -q highmem.q -l h_vmem=20G << EOF
module load MACS2/2.1.0/python.2.7.8
macs2 callpeak -t ${s} -n $out_peaks -g 3e9 
macs2 callpeak -t ${s} -n $out_peaks -g 3e9 --call-summits 
EOF


module load idr/2.0.2/python.3.4.1-atlas-3.11.30 
python3 -c "import idr"
python3 -c "import matplotlib"
idr --samples gm78_rep1_peaks.narrowPeak gm78_rep2_peaks.narrowPeak --input-file-type narrowPeak --output-file gm78_idr --output-file-type narrowPeak --plot 


# and also the IMR90 genotyping
# Getting the data 

module load sra-toolkit 
fastq-dump --outdir fastq --gzip --skip-technical  --readids --dumpbase --split-files --clip SRR_ID

qsub -S /bin/bash -N Download_IMR90 -R y -j y -cwd -q highmem.q -l h_vmem=60G << EOF
module load sra-toolkit
/Users/Andrew/Documents/Greally\ Lab/sratoolkit.2.8.1-mac64/bin/fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip SRR1171946
EOF


fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip SRR1171946

#### Getting the different subunits of NFKb's affinitiy for the locus####
# RelA ChIP-seq rep.1
SRR5579179
# missing the ChIP-seq rep.2 files 
# RelB ChIP-seq rep.1
SRR5579180
# RelB ChIP-seq rep.2
SRR5579181	
# cRel ChIP-seq rep.1
SRR5579182
# cRel ChIP-seq rep.2
SRR5579183
# p50 ChIP-seq rep.1
SRR5579184
# p50 ChIP-seq rep.2
SRR5579185
# p52 ChIP-seq rep.1
SRR5579186
# p52 ChIP-seq rep.2
SRR5579187
 
sra_files=( SRR5579179 SRR5579180 SRR5579181 SRR5579182 SRR5579183 SRR5579184 SRR5579185 SRR5579186 SRR5579187)
printf "%s\n" "${sra_files[@]}"
for sra in "${sra_files[@]}"
do
echo $sra
qsub -S /bin/bash -N Download_NFKB -R y -j y -cwd -l h_vmem=20G -M andrew.johnston@med.einstein.yu.edu -m e << EOF
module load sra-toolkit
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip $sra
EOF
done

zcat SRR5579187_1.fastq.gz  | head -n 1000 | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=1000;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding!";}'
# Phred+33

for sra in "${sra_files[@]}"
do
SampleName="${sra}"
echo $SampleName
qsub -S /bin/bash -N trim_${SampleName} -j y -cwd -l h_vmem=20G << EOF
SampleName="${sra}"
echo $SampleName 
module load trim_galore 
module load fastq_screen/0.4.4/java.1.7.0_67
module unload perl
module load perl/5.16.3/gcc.4.4.7
trim_galore --fastqc ${SampleName}_1.fastq.gz
EOF
done


mv SRR5579179_1_trimmed.fq.gz RelA_1_trimmed.fq.gz 
mv SRR5579180_1_trimmed.fq.gz RelB-rep1_1_trimmed.fq.gz 
mv SRR5579181_1_trimmed.fq.gz RelB-rep2_1_trimmed.fq.gz 
mv SRR5579182_1_trimmed.fq.gz RelC-rep1_1_trimmed.fq.gz 
mv SRR5579183_1_trimmed.fq.gz RelC-rep2_1_trimmed.fq.gz 
mv SRR5579184_1_trimmed.fq.gz p50-rep1_1_trimmed.fq.gz 
mv SRR5579185_1_trimmed.fq.gz p50-rep2_1_trimmed.fq.gz 
mv SRR5579186_1_trimmed.fq.gz p52-rep1_1_trimmed.fq.gz 
mv SRR5579187_1_trimmed.fq.gz p52-rep2_1_trimmed.fq.gz 

####  tst.awk #######
BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}
################################
printf 'TotalReadPair\nDistinctReadPairs\nOneReadPair\nTwoReadPairs\nNRF\nPBC1\nPBC2' > temp.txt


mkdir Bams
mkdir Sams
mkdir Flagstats
mkdir Picard_metrics
module load picard-tools/1.92/java.1.8.0_20 
for file in *trimmed.fq.gz; 
do
SampleName="$(echo ${file} | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N align_$SampleName -R y -j y -cwd -pe smp 20 -R y -l h_vmem=5.6G << EOF
module load bwa/0.7.15/gcc.4.4.7
bwa mem -M -t 20 /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa ${SampleName}_1_trimmed.fq.gz > Sams/${SampleName}.sam
module load picard-tools/1.92/java.1.8.0_20
module load samtools/0.1.19/gcc.4.4.7 

# get unique reads 
samtools view -Sbq 1 -F 4 Sams/${SampleName}.sam | samtools sort -@ 6 - Bams/${SampleName}_UMap
echo "${SampleName}_UMap.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap.bam >> Flagstats/${SampleName}_flagstat.txt
samtools index Bams/${SampleName}_UMap.bam

## mark duplicates 
java -jar $(which MarkDuplicates.jar) INPUT=Bams/${SampleName}_UMap.bam OUTPUT=Bams/${SampleName}_UMap_mkdup.bam METRICS_FILE=Picard_metrics/${SampleName}_UMap_noMT_mkdup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates ASSUME_SORTED=true MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000
echo "${SampleName}_UMap_mkdup.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap_mkdup.bam >> Flagstats/${SampleName}_flagstat.txt

# also remove duplicates with samtools rmdup
samtools rmdup Bams/${SampleName}_UMap_mkdup.bam Bams/${SampleName}_UMap_mkdup_rmdup.bam
echo "${SampleName}_UMap_mkdup_rmdup.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap_mkdup_rmdup.bam >> Flagstats/${SampleName}_flagstat.txt
samtools index Bams/${SampleName}_UMap_mkdup_rmdup.bam

# analyze the NRF and PBC_1 and PBC_2 of the ChIPseq expt
module load bedtools2
bedtools bamtobed -i Bams/${SampleName}_UMap.bam | \
awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,\$6}' | \
grep -v 'chrM' | sort | uniq -c | \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > Flagstats/${SampleName}_PBC_QC_temp.txt

printf 'TotalReadPair\nDistinctReadPairs\nOneReadPair\nTwoReadPairs\nNRF\nPBC1\nPBC2' > temp.txt
awk -f tst.awk Flagstats/${SampleName}_PBC_QC_temp.txt | paste temp.txt - > Flagstats/${SampleName}_PBC_QC.txt 
rm Flagstats/${SampleName}_PBC_QC_temp.txt

# Call peaks 
mkdir Peaks
mkdir Summits
cd Bams
DIR='../Peaks/'
DIR2='../Summits/'
out_peaks="$DIR${SampleName}"
out_summits="$DIR2${SampleName}"
module load MACS2/2.1.0/python.2.7.8
macs2 callpeak -t ${SampleName}_UMap_mkdup_rmdup.bam -n $out_peaks -g 3e9 
macs2 callpeak -t ${SampleName}_UMap_mkdup_rmdup.bam -n $out_summits -g 3e9 --call-summits 
EOF
done


## fixed errors?
# analyze the NRF and PBC_1 and PBC_2 of the ChIPseq expt

for file in *trimmed.fq.gz; 
do
SampleName="$(echo ${file} | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N align_$SampleName -R y -j y -cwd -R y -l h_vmem=30G << EOF
module load bedtools2
bedtools bamtobed -i Bams/${SampleName}_UMap.bam | \
awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,\$6}' | \
grep -v 'chrM' | sort | uniq -c | \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > Flagstats/${SampleName}_PBC_QC_temp.txt

printf 'TotalReadPair\nDistinctReadPairs\nOneReadPair\nTwoReadPairs\nNRF\nPBC1\nPBC2' > temp.txt
awk -f tst.awk Flagstats/${SampleName}_PBC_QC_temp.txt | paste temp.txt - > Flagstats/${SampleName}_PBC_QC.txt 
rm Flagstats/${SampleName}_PBC_QC_temp.txt

# Call peaks 
mkdir Peaks
mkdir Summits
cd Bams
DIR='../Peaks/'
DIR2='../Summits/'
out_peaks="$DIR${SampleName}"
out_summits="$DIR2${SampleName}"
module load MACS2/2.1.0/python.2.7.8
macs2 callpeak -t ${SampleName}_UMap_mkdup_rmdup.bam -n $out_peaks -g 3e9 
macs2 callpeak -t ${SampleName}_UMap_mkdup_rmdup.bam -n $out_summits -g 3e9 --call-summits 
EOF
done

# I don't think that the peak calling worked.

for file in *UMap_mkdup_rmdup.bam; 
do
SampleName="$(echo ${file} | cut -d '_' -f1)"
echo $SampleName
DIR='../Peaks/'
DIR2='../Summits/'
out_peaks="$DIR${SampleName}"
echo $out_peaks
out_summits="$DIR2${SampleName}"
echo $out_summits
qsub -S /bin/bash -N peakCall_$SampleName -R y -j y -q highmem.q -cwd -R y -l h_vmem=15G << EOF
DIR='../Peaks/'
DIR2='../Summits/'
out_peaks="$DIR${SampleName}"
echo $out_peaks
out_summits="$DIR2${SampleName}"
module load MACS2/2.1.0/python.2.7.8
macs2 callpeak -t ${file} -n $out_peaks -g 3e9 
macs2 callpeak -t ${file} -n $out_summits -g 3e9 --call-summits 
EOF
done
## 

## then need to intersect my locus with the various narrowPeak files
vi TBC1D4.bed
##########
chr13	75300370	75300410	TBC1D4	.	.
###########
module load bedtools2

cat myhg19Genes.bed sed 's/^chr//' 

for f1 in *.narrowPeak
do
bedtools intersect -a $f1 -b TBC1D4.bed -nonamecheck
done 