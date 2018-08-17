
## performing the alignment, etc

mkdir Bams
mkdir Sams
mkdir Flagstats
mkdir Picard_metrics

for f1 in *_trimmed.fq.gz
do 
ReadNum="$(echo ${f1} | cut -d '_' -f1 | cut -d '.' -f2)"
echo $ReadNum
SampleName="$(echo ${f1} | cut -d '_' -f1 | cut -d '.' -f1)"
echo $SampleName
if [ $ReadNum != '2' ]
then
Read1="$(echo ${f1})"
echo $Read1	
else
echo ${f1}
qsub -S /bin/bash -N ${SampleName}_align -j y -cwd -pe smp 12 -l h_vmem=5.6G << EOF
module load bwa/0.7.13/gcc.4.4.7
bwa mem -M -t 12 /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa $Read1 ${f1} > Sams/${SampleName}.sam
module load samtools
module load picard-tools/1.92/java.1.8.0_20 

## Select only uniquely mapped reads and sort
samtools view -bq 1 -F 4 Sams/${SampleName}.sam | samtools sort -@ 9 - Bams/${SampleName}_UMap
echo "${f1%.sam}_UMap.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap.bam >> Flagstats/${SampleName}_flagstat.txt
samtools index Bams/${SampleName}_UMap.bam

## Remove Mitochondrial reads 
# Get reads that are not mitochondrial chr
samtools view Bams/${SampleName}_UMap.bam | awk '\$3~"chrM"{print \$1}' | uniq > ${SampleName}.read.list
java -jar $(which FilterSamReads.jar) I=Bams/${SampleName}_UMap.bam O=Bams/${SampleName}_UMap_noMT.bam READ_LIST_FILE=${SampleName}.read.list FILTER=excludeReadList
samtools index Bams/${f1%.sam}_UMap_rmChrM.bam
echo "${SampleName}_UMap_noMT.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap_noMT.bam >> Flagstats/${SampleName}_flagstat.txt
	
## Remove duplicates
java -jar $(which MarkDuplicates.jar) INPUT=Bams/${SampleName}_UMap_noMT.bam OUTPUT=Bams/${SampleName}_UMap_noMT_mkdup.bam METRICS_FILE=Picard_metrics/${SampleName}_UMap_noMT_mkdup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates ASSUME_SORTED=true MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000
echo "${SampleName}_UMap_noMT_mkdup.bam" >> Flagstats/${SampleName}_flagstat.txt
samtools flagstat Bams/${SampleName}_UMap_noMT_mkdup.bam >> Flagstats/${SampleName}_flagstat.txt
samtools index Bams/${SampleName}_UMap_noMT_mkdup.bam
EOF
fi
done

# perform extra samtools rmdup
for f1 in *_UMap_noMT_mkdup.bam;
do
echo $f1
SampleName="$(echo ${f1} | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N ${SampleName}_rmdup -cwd -l h_vmem=20G -j y << EOF
module load samtools/0.1.19/gcc.4.4.7
samtools rmdup ${f1} ${f1%.bam}_rmdup.bam
echo "${SampleName}_UMap_noMT_mkdup_rmdup.bam" >> ../Flagstats/${SampleName}_flagstat.txt
samtools flagstat ${f1%.bam}_rmdup.bam >> ../Flagstats/${SampleName}_flagstat.txt
EOF
done


## Shift the reads for peak calling

mkdir Beds
for s in *_UMap_noMT_mkdup_rmdup.bam;
#for s in gm83_UMap_noMT_mkdup.bam;
do
SampleName="$(echo ${s} | cut -d '_' -f1)"
qsub -S /bin/bash -N ${SampleName}_Shift -cwd -l h_vmem=20G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
module load bedtools2
module load samtools 
samtools view -h -f 0x0040 -b ${s} > ${s%.bam}_read1.bam
bedtools bamtobed -i ${s%.bam}_read1.bam > ../Beds/${SampleName}.bed
awk 'BEGIN {OFS = "\t"} ; {if (\$6 == "+") print \$1, \$2 + 4, \$3 + 4, \$4, \$5, \$6; else print \$1, \$2 - 5, \$3 - 5, \$4, \$5, \$6}' ../Beds/${SampleName}.bed > ../Beds/${SampleName}_shifted.bed
EOF
done

# Rep 1 
mv gm78new* Unused_beds/
mv gm78.bed Unused_beds/
mv gm78_shifted.bed Unused_beds/
mv gm83new* Unused_beds/
mv gm83.bed Unused_beds/
mv gm83_shifted.bed Unused_beds/ 
mv gm87new* Unused_beds/
mv gm87.bed Unused_beds/
mv gm87_shifted.bed Unused_beds/
mv gm92new* Unused_beds/
mv gm92.bed Unused_beds/
mv gm92_shifted.bed Unused_beds/ 
mv gm93new* Unused_beds/
mv gm93.bed Unused_beds/
mv gm93_shifted.bed Unused_beds/


# Rep 2 
mv gm82new* Unused_beds/
mv gm82.bed Unused_beds/
mv gm82_shifted.bed Unused_beds/
mv gm83new* Unused_beds/
mv gm83.bed Unused_beds/
mv gm83_shifted.bed Unused_beds/ 
mv gm92new* Unused_beds/
mv gm92.bed Unused_beds/
mv gm92_shifted.bed Unused_beds/ 
mv gm93new* Unused_beds/
mv gm93.bed Unused_beds/
mv gm93_shifted.bed Unused_beds/


# Clean the bed files
for f1 in *_shifted.bed
do
awk '{if ($1~"chr") print $0; else print "chr"$0}' ${f1} | awk '{if ($2>0) print $0; else print $1"\t"0"\t"$3"\t"$4"\t"$5"\t"$6}' > ${f1%.bed}_chr_noNeg.bed
done

# rep1
cat *_shifted_chr_noNeg.bed > Rep1_shifted.bed
# rep2 
cat *_shifted_chr_noNeg.bed > Rep2_shifted.bed

for s in *shifted_chr_noNeg.bed;
#for s in gm83_UMap_noMT_mkdup.bam;
do
SampleName="$(echo ${s} | cut -d '_' -f1)"
DIR='../Peaks_shift/'
#DIR2='../Peaks_PE_shift/'
out="$DIR${SampleName}"
# outPE="$DIR2${SampleName}_BAMPE"
echo $out
echo $s
qsub -S /bin/bash -N ${SampleName}_Call_Peaks -cwd -l h_vmem=20G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
#Jason's method
macs2 callpeak --nomodel -t ${s} -n $out --nolambda -g 3e9 --keep-dup 'all' --slocal 10000 --call-summits 
EOF
done

for s in Rep2_shifted.bed;
#for s in gm83_UMap_noMT_mkdup.bam;
do
SampleName="$(echo ${s} | cut -d '_' -f1)"
DIR='../Peaks_shift/'
#DIR2='../Peaks_PE_shift/'
out="$DIR${SampleName}"
# outPE="$DIR2${SampleName}_BAMPE"
echo $out
echo $s
qsub -S /bin/bash -N ${SampleName}_Call_Peaks -cwd -q highmem.q -l h_vmem=40G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
#Jason's method
macs2 callpeak --nomodel -t ${s} -n $out --nolambda -g 3e9 --keep-dup 'all' --slocal 10000 --call-summits 
EOF
done

## Get the FRiP
# peaks are written out for each peak
for f1 in *_peaks.narrowPeak;
do
awk '!array[$2,$3]++' $f1 > ${f1%.narrowPeak}_firstPeak.narrowPeak
done

# Clean narrowPeak files
for f1 in *_firstPeak.narrowPeak;
do
awk '{if ($1~"chr") print $0; else print "chr"$0}' ${f1} > ${f1%.narrowPeak}_chr.narrowPeak
done

touch FriP.txt
module load bedtools2
for f1 in *firstPeak_chr.narrowPeak
do
SampleName="$(echo ${f1} | cut -d '_' -f1)"
echo $SampleName >> FriP.txt
echo $f1
bedtools intersect -a ../Beds/${SampleName}_shifted_chr_noNeg.bed -b ${f1} -u | wc -l >> FriP.txt
done

# rep1 
for f1 in Rep1_peaks_firstPeak_chr.narrowPeak
do
SampleName="$(echo ${f1} | cut -d '_' -f1)"
echo $SampleName >> FriP.txt
echo $f1
bedtools intersect -a ../Beds/${SampleName}_shifted.bed -b ${f1} -u | wc -l >> FriP.txt
done

## Get the FRiS 
# clean up the summit files 
for f1 in *_summits.bed;
do
awk '{if ($1~"chr") print $0; else print "chr"$0}' ${f1} > ${f1%.bed}_chr.bed
done

# need to +/- 250 bp to each summit
for f1 in *_summits_chr.bed;
do
awk '{if ($2-250 > 0) print $1"\t"$2-250"\t"$3+250"\t"$4"\t"$5; else print $1"\t0\t"$3+250"\t"$4"\t"$5}' $f1 >${f1%.bed}_250.bed
done 


# What about summits that overlap? 
# for FRiS, we'll ignore
touch FRiS.txt
module load bedtools2
for f1 in *_summits_chr_250.bed
#for f1 in gm77_peaks_chr.narrowPeak;
do
SampleName="$(echo ${f1} | cut -d '_' -f1)"
echo $SampleName >> FRiS.txt
echo $f1
bedtools intersect -a ../Beds/${SampleName}_shifted_chr_noNeg.bed -b ${f1} -u | wc -l >> FRiS.txt
done


# to get IDR summits, it'll be tricky
cp Rep1_summits.bed ../../../Total_peak_analysis/
cp Rep2_summits.bed ../../../Total_peak_analysis/
# in total peak analysis
awk '{print $0"\t."}' Rep1_summits.bed > Rep1_summits_strand.bed
awk '{print $0"\t."}' Rep2_summits.bed > Rep2_summits_strand.bed

# Let's try with summit files
module load idr/2.0.2/python.3.4.1-atlas-3.11.30 
python3 -c "import idr"
python3 -c "import matplotlib"
idr --samples Rep1_summits_strand.bed Rep2_summits_strand.bed --input-file-type bed --output-file replicates_summit_idr --output-file-type bed --plot 

# didn't work.. so doing it like this
# first I'll get IDR peaks because IDR analysis will not work for summits. 
# Need to call peaks without call-summit though

#for s in ../Rep1/Fastq/comb_reads/Beds/Rep1_shifted.bed;
for s in ../Replicates/FastQ/Beds/Rep2_shifted.bed;
do
SampleName="$(echo ${s} | cut -d '/' -f5 | cut -d '_' -f1)"
out="${SampleName}_rmdup"
echo $out
echo $s
done
qsub -S /bin/bash -N ${SampleName}_Call_Peaks -cwd -q highmem.q -l h_vmem=40G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
#Jason's method
macs2 callpeak --nomodel -t ${s} -n $out --nolambda -g 3e9 --keep-dup 'all' --slocal 10000
EOF
done


## How consistent is the newly called non -call-summit peaks with the -call-summit peaks?

bedtools intersect -a Rep1_rmdup_peaks.narrowPeak -b ../Rep1/Fastq/comb_reads/Peaks_shift/Rep1_peaks_firstPeak.narrowPeak -wo | less 
# yes, they're the exact same :) 


#mv _rmdup_peaks.narrowPeak Rep2_rmdup_peaks.narrowPeak
#mv _rmdup_peaks.xls Rep2_rmdup_peaks.xls
#mv _rmdup_summits.bed Rep2_rmdup_summits.bed

module load idr/2.0.2/python.3.4.1-atlas-3.11.30 
python3 -c "import idr"
python3 -c "import matplotlib"
idr --samples Rep1_rmdup_peaks.narrowPeak Rep2_rmdup_peaks.narrowPeak --input-file-type narrowPeak --output-file replicates_summit_idr --output-file-type narrowPeak --plot 

awk '$5 > 540 {print $0}' replicates_summit_idr > replicates_summit_idr_05
wc -l replicates_summit_idr_05
# 108,130 replicates_summit_idr_05

### then look at Get_summit_peak_nonspecific.R to get replicate_summit_peak_nonspecific.bed

# To get IDR summits of +/- 250 bp I will first get the 250 bp summits from the original replicate summit files where overlaps are the highest value 
wc -l Rep2_summits_chr_250.bed
# 342333 Rep2_summits_chr_250.bed

bedtools merge -i Rep2_summits_chr_250.bed -c 4,5 -o collapse >  Rep2_summits_chr_250_merge.bed # 125,706

wc -l Rep2_summits_chr_250_merge.bed
# 286294 Rep2_summits_chr_250_merge.bed

## Do Rep 2 now ### 
wc -l Rep1_summits_chr_250.bed
# 277,900 Rep1_summits_chr_250.bed
bedtools merge -i Rep1_summits_chr_250.bed -c 4,5 -o collapse >  Rep1_summits_chr_250_merge.bed # 125,706

wc -l Rep1_summits_chr_250_merge.bed
# 233,314 Rep1_summits_chr_250_merge.bed

## Now I need to intersect the noOver summits and intersect them with the IDR <.05 peaks 

cp Rep1_summits_chr_250_noOver.bed ../../../../Total_peak_analysis/
cp Rep2_summits_chr_250_noOver.bed ../../../Total_peak_analysis/
############# ########
module load bedtools2
bedtools intersect -a Rep1_summits_chr_250_noOver.bed -b replicate_summit_peak_nonspecific.bed -u -f .5 > Rep1_summits_chr_250_noOver_idr.bed
wc -l Rep1_summits_chr_250_noOver_idr.bed
# 137,148
bedtools intersect -a Rep2_summits_chr_250_noOver.bed -b replicate_summit_peak_nonspecific.bed -u -f 0.5 > Rep2_summits_chr_250_noOver_idr.bed
wc -l Rep2_summits_chr_250_noOver_idr.bed
# 139,442

# add the summit as a column
awk '{print $0"\t"$2+250}' Rep1_summits_chr_250_noOver_idr.bed > Rep1_summits_chr_250_noOver_idr_fix.bed
awk '{print $0"\t"$2+250}' Rep2_summits_chr_250_noOver_idr.bed > Rep2_summits_chr_250_noOver_idr_fix.bed

# Now I need to merge the overlapping IDR summits
cat Rep1_summits_chr_250_noOver_idr_fix.bed Rep2_summits_chr_250_noOver_idr_fix.bed > Summit_idrPeak_chr_250_noOver.bed
bedtools sort -i Summit_idrPeak_chr_250_noOver.bed > Summit_idrPeak_chr_250_noOver_sort.bed
wc -l Summit_idrPeak_chr_250_noOver_sort.bed
# 276,590
bedtools merge -i Summit_idrPeak_chr_250_noOver_sort.bed -c 4,5,6 -o collapse > Summit_idrPeak_chr_250_noOver_sort_merge.txt

## Found the non overlapping merged summits
# accidentally doubled some peaks
wc -l Summit_idrPeak_chr_250_noOver_sort_merge_noOver.txt
# 273842 Summit_idrPeak_chr_250_noOver_sort_merge_noOver.txt
# peaks are written out for each peak

awk '!array[$2,$3]++' Summit_idrPeak_chr_250_noOver_sort_merge_noOver.txt | sort -k 1,1 -k2,2n > Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.txt
wc -l Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.txt
# 145,780

awk '{print $1"\t"$2"\t"$3"\tSummit_"NR"\t"$5"\t."}' Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.txt > Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed

# Had to remove .5s from the Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed in R
#  grep "2300853.5" Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed


bedtools intersect -a Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed -b replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed -u | wc -l
# 139,588

bedtools intersect -a Summit_idrPeak_chr_250_noOver_sort_merge_noOver_fix.bed -b replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed -u > Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack.bed

# using  Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed

# Now I can quantify these peaks!
# Rep1 
for f1 in ~/17_member/ATACseq/Rep1/Fastq/comb_reads/Beds/*_shifted_chr_noNeg.bed 
#for f1 in ~/17_member/ATACseq/Rep1/Fastq/comb_reads/Beds/gm77_shifted_chr_noNeg.bed 
do
SampleName="$(echo ${f1} | cut -d '/' -f12 | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N Quant_summit_${SampleName} -cwd -l h_vmem=10G -j y << EOF
module load bedtools2
bedtools intersect -a ${f1} -b ~/17_member/ATACseq/Total_peak_analysis/Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -wo | awk '{print \$10}' | cut -d '_' -f2 | sort -n | uniq -c > ${SampleName}_rep1_readPeaks
EOF
done

for f1 in ~/17_member/ATACseq/Replicates/FastQ/Beds/*_shifted_chr_noNeg.bed 
#for f1 in ~/17_member/ATACseq/Replicates/FastQ/Beds/gm77_shifted_chr_noNeg.bed 
do
SampleName="$(echo ${f1} | cut -d '/' -f11 | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N Quant_summit_${SampleName} -cwd -l h_vmem=10G -j y << EOF
module load bedtools2
bedtools intersect -a ${f1} -b ~/17_member/ATACseq/Total_peak_analysis/Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -wo | awk '{print \$10}' | cut -d '_' -f2 | sort -n | uniq -c > ${SampleName}_rep2_readPeaks
EOF
done

# changing the names of some of the files for eaiser manipulation 
mv gm82comb_rep2_readPeaks gm82_rep2_readPeaks
mv gm83comb_rep2_readPeaks gm83_rep2_readPeaks
mv gm92comb_rep2_readPeaks gm92_rep2_readPeaks
mv gm93comb_rep2_readPeaks gm93_rep2_readPeaks

mv gm78comb_rep1_readPeaks gm78_rep1_readPeaks
mv gm83comb_rep1_readPeaks gm83_rep1_readPeaks
mv gm92comb_rep1_readPeaks gm92_rep1_readPeaks
mv gm93comb_rep1_readPeaks gm93_rep1_readPeaks
mv gm87comb_rep1_readPeaks gm87_rep1_readPeaks

grep -e "e" Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack.bed

# For the CQN normalization, I need to have the GC content and length for each gene.
module load bedtools2
bedtools nuc -fi /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa -bed Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed > Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_GCContent.txt

awk '{print $4"\t"$8"\t"$15}' Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_GCContent.txt > Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_GCContent_DF.txt

awk '{print $3}' Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_GCContent_DF.txt | sort  | uniq -c 
# all of the summits used are 500 bp long

#Performing the Quantification of the summits in R in ~/17_member/ATACseq/Quant_summit_2



###############################################################
######                                                  #######
######                 Getting individual IDR peak      #######
###############################################################

# To get the peaks that at least 3 kids have to filter tested summits
# I need to first get the peaks that the kids call 


# getting the peaks for each sample for each replicate
for s in ../Beds/*shifted_chr_noNeg.bed;
do
SampleName="$(echo ${s} | cut -d '/' -f3 | cut -d '_' -f1)"
out="${SampleName}"
echo $out
echo $s
qsub -S /bin/bash -N ${SampleName}_Call_Peaks -cwd -q highmem.q -l h_vmem=15G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
#Jason's method
macs2 callpeak --nomodel -t ${s} -n $out --nolambda -g 3e9 --keep-dup 'all' --slocal 10000
EOF
done



# Now to perform the IDR analysis 
cd ~/17_member/ATACseq/IDR_peaks/IDR_2017/

module load idr/2.0.2/python.3.4.1-atlas-3.11.30 
python3 -c "import idr"
python3 -c "import matplotlib"
idr --samples Rep1_rmdup_peaks.narrowPeak Rep2_rmdup_peaks.narrowPeak --input-file-type narrowPeak --output-file replicates_summit_idr --output-file-type narrowPeak --plot 

# changed rep2 so that code would work in ~/17_member/ATACseq/Replicates/FastQ/Peaks_Final_rmdup
mv gm82comb_peaks.narrowPeak gm82_peaks.narrowPeak
mv gm87_peaks.narrowPeak gm87comb_peaks.narrowPeak
mv gm78_peaks.narrowPeak gm78comb_peaks.narrowPeak

for f1 in ../../Rep1/Fastq/comb_reads/Peaks_Final_rmdup/*narrowPeak
for f1 in ../../Rep1/Fastq/comb_reads/Peaks_Final_rmdup/gm82_peaks.narrowPeak ../../Rep1/Fastq/comb_reads/Peaks_Final_rmdup/gm87comb_peaks.narrowPeak ../../Rep1/Fastq/comb_reads/Peaks_Final_rmdup/gm78comb_peaks.narrowPeak;
do
SampleName="$(echo ${f1} | cut -d '/' -f7 | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N ${SampleName}_idr -cwd -l h_vmem=20G -j y << EOF
module load idr/2.0.2/python.3.4.1-atlas-3.11.30 
python3 -c "import idr"
python3 -c "import matplotlib"
idr --samples ${f1} ../../Replicates/FastQ/Peaks_Final_rmdup/${SampleName}_peaks.narrowPeak --input-file-type narrowPeak --output-file ${SampleName}_peaks_idr --output-file-type narrowPeak --plot 
EOF
done

# get number of peaks that are IDR < 0.05 
for f1 in *_peaks_idr
do
awk '$5 > 540 {print $0}' $f1 | wc -l
done

# Make < 0.05 peak list 
for f1 in *_peaks_idr
do
awk '$5 > 540 {print $0}' $f1 > ${f1}_05
done

# how many of them intersect replicate_summit_peak_nonspecific.bed
# first clean up the replicate summit peak
sort -k 1,1 -k2,2n replicate_summit_peak_nonspecific.bed | awk '{print $1"\t"$2"\t"$3"\tPeak_"NR"\t"$5"\t."}' > replicate_summit_peak_nonspecific_sort.bed 

# Cleaning up more in the R file to get canonical chr
wc -l replicate_summit_peak_nonspecific_sort_Canon.bed 
# 103691 replicate_summit_peak_nonspecific_sort_Canon.bed

module load bedtools2 
bedtools intersect -a replicate_summit_peak_nonspecific_sort_Canon.bed -b /home/greally-lab/indexes/consensusBlacklist_hg38.bed -v > replicate_summit_peak_nonspecific_sort_Canon_FilMappBlack.bed

wc -l replicate_summit_peak_nonspecific_sort_Canon_FilMappBlack.bed 
#103481 replicate_summit_peak_nonspecific_sort_Canon_FilMappBlack.bed

module load bedtools2 
bedtools intersect -a replicate_summit_peak_nonspecific_sort_Canon_FilMappBlack.bed -b Mito_peak_blacklist.bed -v > replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed

wc -l replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed
# 103,476 replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed


for f1 in *_peaks_idr_05;
do
SAMPLE="$(echo ${f1} | cut -d '_' -f1)"
echo $SAMPLE
bedtools intersect -a ../../Total_peak_analysis/replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed -b ${f1} -wo | awk '{print $4}' | sort | uniq | awk '{print $1"\t1"}' > ${SAMPLE}_inter_master.txt
done

mv gm78comb_inter_master.txt gm78_inter_master.txt
mv gm83comb_inter_master.txt gm83_inter_master.txt
mv gm87comb_inter_master.txt gm87_inter_master.txt
mv gm92comb_inter_master.txt gm92_inter_master.txt
mv gm93comb_inter_master.txt gm93_inter_master.txt

awk '{print $4}' ../../Total_peak_analysis/replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed > Master_peak_names.txt
wc -l Master_peak_names.txt
# 103476 Master_peak_list.txt

wc -l Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack.bed 
# 139588
bedtools merge -i Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack.bed | wc -l
# 139,581

bedtools merge -i Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack.bed -c 4,5 -o collapse >  Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_merge.bed 
# made Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix.bed
# in R 

bedtools sort -i Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix.bed > Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed
wc -l Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed 
# 139588 Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed

bedtools merge -i Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed | wc -l
# 139,588

wc -l replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed 
# 103,476 replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed

bedtools merge -i replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed | wc -l # 103474

# Use this to do the quantification. 

# need to figure out which peaks correspond to which summits
less Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed

module load bedtools2
bedtools intersect -a Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -b replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed -wo > Summit_Summit_peak_inter.txt
wc -l Summit_Summit_peak_inter.txt 
# 139,700 Summit_Summit_peak_inter.txt
less Summit_Summit_peak_inter.txt

bedtools intersect -a Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -b replicate_summit_peak_nonspecific_sort_Canon_FilMappMitoBlack.bed -wo | awk '{print $4}' | sort | uniq -d | wc -l

## to do the chrQTL analysis, I need to intersect the summits with the haplotype information
# pwd = ~/17_member/QTLs_2017/ATAC_summit_QTL_2
module load bedtools2
bedtools intersect -a ../../ATACseq/Total_peak_analysis/Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -b ../haplotype_grch38_named_4intersect.txt -wo > ATAC_summit_Plat_haplo_interesect.txt

wc -l ATAC_summit_Plat_haplo_interesect.txt
#139,301 ATAC_summit_Plat_haplo_interesect.txt
# only 287 summits did not overlap a haplotype

bedtools intersect -a ../../ATACseq/Total_peak_analysis/Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -b ../haplotype_grch38_named_4intersect.txt -wa | uniq -d
# no summits overlap multiple haplotypes





