
Can I estimate the copy number of EBV from the methylation data? 

For all of the samples, the `BS-TAG-Sample-77.filterFillIn.bam` file is the final bam that should be used. does it have EBV ?
```bash
cd /gs/gsfs0/users/anjohnst/17_member/BS-seq/gm77/analysis_hg38
module load samtools
samtools view  BS-TAG-Sample-77.filterFillIn.bam | awk '$3=="chrEBV"' | less 

tail gm77_idxstats.txt 
```
There are reads aligning to the EBV genome. So now we need to think about how to calculate the average coverage. Additionally, it should be noted that the EBV reference genome used is for the B95-8 strain, which has a 12kb deletion when comapred to wildtype EBV (Raji strain). Because we did not align with the addition of this 12kb region, we cannot tell which LCLs within the family have the additional deletion which would provide evidence toward the existence of other EBV strain infections within the cell lines. 

from Santpere, G., Darre, F., Blanco, S., Alcami, A., Villoslada, P., Mar Albà, M., & Navarro, A. (2014). Genome-wide analysis of wild-type Epstein-Barr virus genomes derived from healthy individuals of the 1,000 Genomes Project. Genome Biology and Evolution, 6(4), 846–860. http://doi.org/10.1093/gbe/evu054
```
" However, B95-8 presents a specific 12-kb deletion in its genome that was spontaneously produced in the laboratory (Skare et al. 1982), allowing us to distinguish between natural and artificial EBV strains."

"a composite of the B95-8 strain and 12 kb from the Raji strain to correct the nonnatural B95-8-specific deletion from position 139,724 to 151,554."

"Previous to mapping, we masked the reference EBV genome using RepeatMasker (Smit et al. 1996–2010) to remove low-complexity regions. We also masked a large region of repeats between positions 12,001 and 35,355." 
```

module load R
java -jar picard.jar CollectWgsMetricsWithNonZeroCoverage \
       I=input.bam \
       O=collect_wgs_metrics.txt \
       CHART=collect_wgs_metrics.pdf  \
       R=reference_sequence.fasta 


samtools view  BS-TAG-Sample-77.filterFillIn.bam | tail


I decided not to estimate the CN from the WGS from illumina/1000G because the EBV CN changes over passages. 

I could try to estimate from the X-WGBS data.. but I haven't 

I have the 171823 basepairs in my assembly which includes the raji deletion (https://www.ncbi.nlm.nih.gov/nuccore/NC_007605.1?report=genbank). Let's isolate the EBC chromosome from the bam files 

The deletion is `chrEBV	139724	151554	+	EBV_del`

```bash
module load RepeatMasker/4.0.7/perl.5.22.1



```



I made a bed file of just the deletion and one chrEBV coords without the deletion
deletion (chrEBV:139724-151554)
```
chrEBV	139724	151554	+	EBV_del
```
chrEBV minus the del
```
chrEBV  1       139724  +       EBV_1
chrEBV  151554  171823  +       EBV_2
```

Repeat regions:
```
chrEBV	7368	8080
chrEBV:12,173-35,276
chrEBV:38,260-39,762
chrEBV:40,499-41,698
chrEBV:95,776-96,909
chrEBV:169,472-172,117
```

chrEBV minus the del and repeats: EBV_interval_noRep.bed
```
chrEBV  1       7368	+	EBV_1
chrEBV	8080	12173	+	EBV_2
chrEBV	35276	38260	+	EBV_3
chrEBV	39762	40499	+	EBV_4
chrEBV	41698	95776	+	EBV_5
chrEBV	96909	139724	+	EBV_6
chrEBV	151554  169472	+	EBV_7
```



I made picard interval files

```bash
module load picard/2.17.1/java.1.8.0_20 
java -jar $(which picard.jar) BedToIntervalList \
      I=EBV_interval.bed \
      O=EBV.interval_list \
      SD=BStag_77_filterFillIn_EBV.bam

java -jar $(which picard.jar) BedToIntervalList \
      I=EBV_del_interval.bed \
      O=EBV_del.interval_list \
      SD=BStag_77_filterFillIn_EBV.bam

java -jar $(which picard.jar) BedToIntervalList \
      I=EBV_interval_noRep.bed \
      O=EBV_interval_noRep.interval_list \
      SD=BStag_77_filterFillIn_EBV.bam

```

Next I wanted to only get depth of the canonical chromosomes

Created a file `Canon_coords.bed`
```
chr1	1	248956422
chr2 	1	242193529
chr3 	1	198295559
chr4 	1	190214555
chr5	1	181538259
chr6	1	170805979
chr7	1	159345973
chr8	1	145138636
chr9	1	138394717
chr10	1	133797422
chr11	1	135086622
chr12	1	133275309
chr13	1	114364328
chr14	1	107043718
chr15	1	101991189
chr16	1	90338345
chr17	1	83257441
chr18	1	80373285
chr19	1	58617616
chr20	1	64444167
chr21	1	46709983
chr22	1	50818468
```
Now to subsract out blacklisted regions 
```bash
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gunzip -d hg38.blacklist.bed.gz

module load bedtools2
bedtools subtract -a Canon_coords.bed -b hg38.blacklist.bed > Canon_coords_masked.bed
```
I then added the + and unique interval name and made file `Canon_coords_masked_fix.bed`
```bash
awk '{OFS="\t"; print $0,"+","region_"NR}' Canon_coords_masked.bed > Canon_coords_masked_fix2.bed
```

Then I needed to make an interval for picard
```bash
java -jar $(which picard.jar) BedToIntervalList \
      I=Canon_coords_masked_fix2.bed \
      O=Canon_coords_masked.interval_list \
      SD=BStag_77_filterFillIn_EBV.bam
```


java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=BStag_93_filterFillIn_EBV.bam O=test_EBV_wgs_metrics.txt CHART=test_EBV_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV.interval_list COVERAGE_CAP=999999999 LOCUS_ACCUMULATION_CAP=999999999

```bash
for f1 in gm78 gm79 gm80 gm81 gm82 gm83 gm84 gm85 gm86 gm87 gm88 gm89 gm90 gm91 gm92 gm93;
module load picard/2.17.1/java.1.8.0_172

for f1 in gm77;
for f1 in gm78 gm79 gm80 gm81 gm82 gm83 gm84 gm85 gm86 gm87 gm88 gm89 gm90 gm91 gm92 gm93;
for f1 in gm*;
do
SampleName="$(echo ${f1} | cut -c 3-4)"
echo $SampleName
dir="$(echo ${f1}"/analysis_hg38/")"
echo $dir
if [ $SampleName != '84' ] && [ $SampleName != '89' ]
then
echo "Not 84 or 89"
qsub -S /bin/bash -N getEBV_${SampleName} -q highmem.q -l h_vmem=300G -j y -cwd << EOF
module load samtools/1.2/gcc.4.4.7
module load picard/2.17.1/java.1.8.0_20
module load R/3.4.0

samtools view -b ${dir}BS-TAG-Sample-${SampleName}.filterFillIn.bam chrEBV > EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam
samtools index EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam

java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam O=EBV_Bams/${SampleName}_EBV_wgs_metrics.txt CHART=EBV_Bams/${SampleName}_EBV_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_Bams/EBV.interval_list COVERAGE_CAP=999999999 LOCUS_ACCUMULATION_CAP=999999999

java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam O=EBV_Bams/${SampleName}_EBVdel_wgs_metrics.txt CHART=EBV_Bams/${SampleName}_EBVdel_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_Bams/EBV_del.interval_list COVERAGE_CAP=999999999 LOCUS_ACCUMULATION_CAP=999999999

java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=${dir}BS-TAG-Sample-${SampleName}.filterFillIn.bam O=EBV_Bams/${SampleName}_canon_wgs_metrics.txt CHART=EBV_Bams/${SampleName}_canon_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_Bams/Canon_coords_masked.interval_list
EOF
else
echo "Yes, 84 or 89"
qsub -S /bin/bash -N getEBV_${SampleName} -q highmem.q -l h_vmem=300G -j y -cwd << EOF
module load samtools/1.2/gcc.4.4.7
module load picard/2.17.1/java.1.8.0_20
module load R/3.4.0

samtools view -b ${dir}BSTag-Sample-${SampleName}-reprep.filterFillIn.bam chrEBV > EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam
samtools index EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam

java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam O=EBV_Bams/${SampleName}_EBV_wgs_metrics.txt CHART=EBV_Bams/${SampleName}_EBV_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_Bams/EBV.interval_list COVERAGE_CAP=999999999 LOCUS_ACCUMULATION_CAP=999999999

java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=EBV_Bams/BStag_${SampleName}_filterFillIn_EBV.bam O=EBV_Bams/${SampleName}_EBVdel_wgs_metrics.txt CHART=EBV_Bams/${SampleName}_EBVdel_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_Bams/EBV_del.interval_list COVERAGE_CAP=999999999 LOCUS_ACCUMULATION_CAP=999999999

java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=${dir}BSTag-Sample-${SampleName}-reprep.filterFillIn.bam O=EBV_Bams/${SampleName}_canon_wgs_metrics.txt CHART=EBV_Bams/${SampleName}_canon_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_Bams/Canon_coords_masked.interval_list
EOF
fi
done

module load R/3.4.0/gcc.4.7.4
```

Analyzing the CN in R
```R
all_cov_nonZero <-data.frame(row.names=c("WG","EBV","EBV_del"))
all_cov <-data.frame(row.names=c("WG","EBV","EBV_del"))
for (sample in seq(77,93,1)){
cov_wg <- read.table(paste(sample, "_canon_wgs_metrics.txt", sep=""), skip=6, nrow=2, header=T)
cov_ebv_del <- read.table(paste(sample,"_EBVdel_wgs_metrics.txt",sep=""), skip=6, nrow=2, header=T)
cov_ebv <- read.table(paste(sample,"_EBV_wgs_metrics.txt",sep=""), skip=6, nrow=2, header=T)

cov_nonZero <- rbind(cov_wg[2,],cov_ebv[2,],cov_ebv_del[2,])
cov_all <- rbind(cov_wg[1,],cov_ebv[1,],cov_ebv_del[1,])
 # row.names(cov_nonZero) <- c(paste(sample,"_WG", sep=""),paste(sample,"_EBV", sep=""))
df <- data.frame(cov_nonZero[,3])
colnames(df) <- sample
all_cov_nonZero <- cbind(all_cov_nonZero, df)

df2 <- data.frame(cov_all[,3])
colnames(df2) <- sample
all_cov <- cbind(all_cov, df2)

if (cov_ebv_del[1,3]>1) {
	print(paste(sample, "has high coverage over EBV deletion"))
}
}

all_cov_nonZero
all_cov
CN <- all_cov_nonZero[2,]/all_cov_nonZero[1,]
```

78=25CN
80=21.7CN

24.600713=2 copies	245.417036 	 8.341176
```
all_cov[3,]
              77       78       79       80       81       82       83       84
EBV_del 0.179797 0.135503 0.100169 0.186306 0.109721 0.142096 0.220541 0.056889
              85       86       87      88       89       90       91       92
EBV_del 0.110059 0.189096 0.224345 0.13153 0.055959 0.143787 0.271851 0.331868
              93
EBV_del 0.173457
```

I need to downsample the BAMs so that I can accurately estimate the coverage

module load sambamba/0.5.4


for f1 in BStag_*_filterFillIn_EBV.bam
do
SampleName="$(echo ${f1} | cut -d '_' -f2)"
echo $SampleName 
qsub -S /bin/bash -N sub_${SampleName} -l h_vmem=2G -j y -cwd << EOF
module load samtools/1.2/gcc.4.4.7
samtools view -b -s 0.004 $f1 > ${f1%.bam}_sub4.bam
EOF
done

for f1 in BStag_*_filterFillIn_EBV.bam
do
SampleName="$(echo ${f1} | cut -d '_' -f2)"
echo $SampleName 
qsub -S /bin/bash -N getEBV_${SampleName} -q highmem.q -l h_vmem=950G -j y -cwd << EOF
module load picard/2.17.1/java.1.8.0_20
module load R/3.4.0
java -jar $(which picard.jar) CollectWgsMetricsWithNonZeroCoverage I=BStag_${SampleName}_filterFillIn_EBV.bam O=${SampleName}_EBV_noCAP_wgs_metrics.txt CHART=${SampleName}_EBV_noCAP_wgs_metrics.pdf R=/home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa INTERVALS=EBV_interval_noRep.interval_list COVERAGE_CAP=5000 LOCUS_ACCUMULATION_CAP=5000
EOF
done

qrsh -l h_vmem=1300G -q highmem.q


```R
all_cov_nonZero <-data.frame(row.names=c("WG","EBV"))
all_cov <-data.frame(row.names=c("WG","EBV"))
for (sample in seq(77,93,1)){
cov_wg <- read.table(paste(sample, "_canon_wgs_metrics.txt", sep=""), skip=6, nrow=2, header=T)
cov_ebv <- read.table(paste(sample,"_EBV_noCAP_wgs_metrics.txt",sep=""), skip=6, nrow=2, header=T)

cov_nonZero <- rbind(cov_wg[2,],cov_ebv[2,])
cov_all <- rbind(cov_wg[1,],cov_ebv[1,])
 # row.names(cov_nonZero) <- c(paste(sample,"_WG", sep=""),paste(sample,"_EBV", sep=""))
df <- data.frame(cov_nonZero[,3])
colnames(df) <- sample
all_cov_nonZero <- cbind(all_cov_nonZero, df)

df2 <- data.frame(cov_all[,3])
colnames(df2) <- sample
all_cov <- cbind(all_cov, df2)
}

all_cov_nonZero
all_cov
CN <- all_cov_nonZero[2,]/all_cov_nonZero[1,]

write.table(CN, "Estimated_EBV_CN.txt", sep = "\t", append = F, row.names = T, col.names = T, quote = F)

```

78=25CN
80=21.7CN

24.600713=2 copies	245.417036 	 8.341176
```
all_cov[3,]
              77       78       79       80       81       82       83       84
EBV_del 0.179797 0.135503 0.100169 0.186306 0.109721 0.142096 0.220541 0.056889
              85       86       87      88       89       90       91       92
EBV_del 0.110059 0.189096 0.224345 0.13153 0.055959 0.143787 0.271851 0.331868
              93
EBV_del 0.173457
```

