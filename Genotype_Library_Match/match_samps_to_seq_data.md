# Verifying family sequencing data matches.
## Andrew D Johnston
### 03/03/2018

Checking that the genotype matches between the family and their respective sequencing data. 

the vcf file will be the platinum genomes vcf found here: `/gs/gsfs0/users/anjohnst/17_member/Platinum_Genomes/files/GenotypeFiles/phg000812.v1.PlatinumGenomes.genotype-calls-vcf.c1/High_Confidence_Calls_HG38.vcf`

First let's prep the vcf file (it's contigs are without the "chr" prefix)

```bash
mkdir 17_member/Match_files
cd 17_member/Match_files
mkdir Genotype
cd Genotype
cp ~/17_member/Platinum_Genomes/files/GenotypeFiles/phg000812.v1.PlatinumGenomes.genotype-calls-vcf.c1/High_Confidence_Calls_HG38.vcf .

module load htslib/1.2.1/gcc.4.4.7 
bgzip High_Confidence_Calls_HG38.vcf
tabix -p vcf High_Confidence_Calls_HG38.vcf.gz 
 # check to see if the tbi is readable by QTLtools
module load bcftools/1.3.1/gcc.4.4.7
bcftools view High_Confidence_Calls_HG38.vcf.gz | less -S
```

Now let's begin matching the various sequencing samples that we have. ATACseq rep1, ATACseq rep2, RNAseq, sRNAseq, and WGBS (the WGBS matching may not be ideal since it is bisulphite converted).

### ATACseq matching

Starting with ATACseq datasets. The first replicate is found in `/gs/gsfs0/users/anjohnst/17_member/ATACseq/Rep1/Fastq/comb_reads/Bams`


```bash
cd ~/17_member/ATACseq/Rep1/Fastq/comb_reads/Bams

for f1 in *_mkdup.bam;
do
SAMPLE="$(echo ${f1} | cut -d '_' -f1)"
echo $SAMPLE
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  ~/17_member/Match_files/Genotype/High_Confidence_Calls_HG38.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ~/17_member/Match_files/ATACR1_match/${SAMPLE}.match.txt
EOF
done

cd ~/17_member/Match_files/ATACR1_match/

for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -d '_' -f1)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript ../plot_sample_match.R $f1
EOF
done
```

Now for the second replicate which is found in `~/17_member/ATACseq/Replicates/FastQ/Bams`

```bash
cd ~/17_member/ATACseq/Replicates/FastQ/Bams

for f1 in *_mkdup.bam;
do
SAMPLE="$(echo ${f1} | cut -d '_' -f1)"
echo $SAMPLE
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  ~/17_member/Match_files/Genotype/High_Confidence_Calls_HG38.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ~/17_member/Match_files/ATACR2_match/${SAMPLE}.match.txt
EOF
done

cd ~/17_member/Match_files/ATACR2_match/
cp ../ATACR1_match/plot_sample_match.R

for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
echo $SAMPLE
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done
```

### RNAseq Matching 


The mapped RNAseq files can be found in `~/17_member/RNA-seq/Mapped_GenCode_EBV_decoy_1pass`

```bash
mkdir ../RNAseq_match
cd ~/17_member/RNA-seq/Mapped_GenCode_EBV_decoy_1pass

for f1 in *Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=12G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  ~/17_member/Match_files/Genotype/High_Confidence_Calls_HG38.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ~/17_member/Match_files/RNAseq_match/${SAMPLE}.match.txt
EOF
done

cd ~/17_member/Match_files/RNAseq_match/
cp ../ATACR2_match/plot_sample_match.R .

for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done
```

### sRNAseq Matching 


The mapped RNAseq files can be found in `~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy`

```bash
mkdir ../sRNAseq_match
cd ~/17_member/miRNA/Run2/Combined_reads/Mapped_GenCode_EBV_decoy

for f1 in *Aligned.sortedByCoord.out.filtered.primary.bam;
do
SAMPLE="$(echo ${f1} | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=12G -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  ~/17_member/Match_files/Genotype/High_Confidence_Calls_HG38.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ~/17_member/Match_files/sRNAseq_match/${SAMPLE}.match.txt
EOF
done

cd ~/17_member/Match_files/sRNAseq_match/
cp ../ATACR2_match/plot_sample_match.R .

for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done
```

### WGBS match

I am somewhat skeptical that this will work because of the bisulphite conversion. However, I am hoping that the the correct sample will have a higher proportion of correct homo/hetero variants when compared to the others. 

These bam file are all in their sample folders within `~/17_member/BS-seq/`. 

```bash
mkdir ../WGBS_match
cd ~/17_member/BS-seq/
ls gm*/analysis_hg38/*.filterFillIn.bam 

for f1 in gm*/analysis_hg38/*.filterFillIn.bam ;
do
SAMPLE="$(echo ${f1} | cut -d '/' -f3 | cut -d '-' -f4 | cut -d '.' -f1)"
echo $SAMPLE
if [ $SAMPLE != 'reprep' ]
then
echo $f1
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=100G -q highmem.q -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  ~/17_member/Match_files/Genotype/High_Confidence_Calls_HG38.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ~/17_member/Match_files/WGBS_match/${SAMPLE}.match.txt
EOF
else 
SAMPLE="$(echo ${f1} | cut -d '/' -f3 | cut -d '-' -f3)"
echo $SAMPLE
echo $f1
qsub -S /bin/bash -N match_${SAMPLE} -cwd -l h_vmem=100G -q highmem.q -j y << EOF
module load QTLtools/1.1/gcc.4.7.4
QTLtools mbv --bam $f1 --vcf  ~/17_member/Match_files/Genotype/High_Confidence_Calls_HG38.vcf.gz --filter-mapping-quality 0 --filter-keep-duplicates --filter-base-quality 5 --filter-minimal-coverage 10 --out ~/17_member/Match_files/WGBS_match/${SAMPLE}.match.txt
EOF
fi
done

cd ~/17_member/Match_files/WGBS_match/
cp ../sRNAseq_match/plot_sample_match.R .

for f1 in *.match.txt;
do
SAMPLE="$(echo ${f1} | cut -d '.' -f1)"
qsub -S /bin/bash -N plotMatch_${SAMPLE} -cwd -l h_vmem=1G -j y << EOF
module load R/3.4.0/gcc.4.7.4
Rscript plot_sample_match.R $f1
EOF
done
```





