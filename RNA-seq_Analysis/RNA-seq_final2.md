# Final quantification and eQTL calling for the 17 member family 
## Andrew D. Johnston
### 7/1/2018

Here is a detailed workflow of how I generated the eQTLs


1. Building the Kallisto genome indexes
	
I have the current one (all transcripts+EBV), but I want to create two more: protein coding+EBV and protein coding alone. 

I use CGAT tools to filter the gtf files. Here is the download and install information
```bash
 # download the tools
curl -O https://raw.githubusercontent.com/CGATOxford/cgat/master/install-CGAT-tools.sh

 # install the tools
./install-CGAT-tools.sh --production 

source /gs/gsfs0/users/anjohnst/cgat-install/conda-install/bin/activate cgat-s
```

Now to filter for just protein coding and linc "genes"
```bash
 # filter for prot coding genes  
cgat gtf2gtf --method=filter --filter-method=proteincoding -I /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/gencode.v24.primary_assembly.annotation.EBV.gtf > gencode.v24.primary_assembly.annotation.prot.gtf

 # filter for lincRNAs
cgat gtf2gtf --method=filter --filter-method=lincrna -I /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/gencode.v24.primary_assembly.annotation.EBV.gtf > gencode.v24.primary_assembly.annotation.lincRNA.gtf

 # Combine
cat gencode.v24.primary_assembly.annotation.prot.gtf gencode.v24.primary_assembly.annotation.lincRNA.gtf > gencode.v24.primary_assembly.annotation.prot_lincRNA.gtf

 # check
awk '$3=="gene" {print $12}' gencode.v24.primary_assembly.annotation.prot_lincRNA.gtf | sort | uniq -c 
 #   7674 "lincRNA";
 #  19844 "protein_coding";

 # parse out the EBV genes (there are 91)
grep "EBV" /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/gencode.v24.primary_assembly.annotation.EBV.gtf > gencode.v24.primary_assembly.annotation.EBV_only.gtf

 # combine the EBV with the protein coding and lincRNA
cat gencode.v24.primary_assembly.annotation.prot_lincRNA.gtf gencode.v24.primary_assembly.annotation.EBV_only.gtf > gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV.gtf

 # remove chrM and HLA genes and 1 transcript duplicate (due to linc and protein coding filter)
grep -v "chrM" gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV.gtf >  gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM.gtf

grep -v "HLA-" gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM.gtf > gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA.gtf

 # find the transcripts that were duplicated
awk '$3=="transcript" {print $12}' gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA.gtf | sort | uniq -c | awk '$1>=2' | less

grep -Ev "ENST00000633294.1|ENST00000306731.4|ENST00000623205.2|ENST00000623295.1" gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA.gtf > gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix.gtf 

grep -E "ENST00000633294.1|ENST00000306731.4|ENST00000623205.2|ENST00000623295.1" gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA.gtf | awk 'NR>12' > ENST_dup.gtf 

cat gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix.gtf ENST_dup.gtf > gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_presort.gtf

 # sort the gtf
cgat gtf2gtf --method=sort --sort-order=gene  -I gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_presort.gtf > gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_sort.gtf

 # create a fasta file for Kallisto
module load cufflinks/2.2.1/gcc.4.4.7
gffread gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_sort.gtf -g /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fa -w gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_sort.fasta​
```

Building Kallisto index

```bash
 # Build Indexes
qsub -S /bin/bash -N Build_kallisto_protlincEBV -cwd -l h_vmem=20G -j y << EOF
module load kallisto
kallisto index -i gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_sort.idx gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_sort.fasta​
EOF
```
v0.42.5
Running kallisto
```bash
 # Alignment and quantification of transcriptomes
for f1 in  ../Trimmed_Reads/Cat_Fq/*.fq.gz
do 
FILE="$(echo ${f1} | cut -d '/' -f4)"
ReadNum="$(echo ${FILE} | cut -d '.' -f2)"
echo $ReadNum
SampleName="$(echo ${FILE} | cut -d '.' -f1)"
echo $SampleName
if [ $ReadNum != '2' ]
then
Read1="$(echo ${f1})"
echo $Read1	
else
echo ${f1}
mkdir ${SampleName}_kallisto_protlincEBV
qsub -S /bin/bash -N ${SampleName}_quant_protlincEBV -j y -cwd -pe smp 10 -l h_vmem=5.6G << EOF
module load kallisto
kallisto quant -i gencode.v24.primary_assembly.annotation.prot_lincRNA_EBV_noM_noHLA_fix_sort.idx --bias -o ${SampleName}_kallisto_protlincEBV -b 100 -t 10 $Read1 ${f1}  
EOF
fi
done
```
--bias -b 100
2. Sleuth to gene eQTL Calling


3. Inegrating with other information

4. Curating final SNP list with rsIDs
downloaded most recent SNP table from UCSC 


awk '{OFS="\t"; print $2,$3,$4,$5,$7,$9,$11,$12}' SNP150.txt > SNP150_table.txt
awk 'NR>1 {OFS="\t"; print $2,$3-1,$4+1,$5}' SNP150.txt > SNP150.bed
sort -k 1,1 -k2,2n eQTL_peakIM_check_vars2.txt > eQTL_peakIM_check_vars2.bed
sort -k 1,1 -k2,2n SNP150.bed > SNP150_sort.bed
grep "rs61960554" SNP150_sort.bed 

awk '{print $1}' SNP150_sort.bed | uniq
module load bedtools2
bedtools intersect -a  eQTL_peakIM_check_vars2.bed -b SNP150_sort.bed -wo | less 
bedtools intersect -a  eQTL_sumIM_check_vars2.txt -b SNP150.bed -wo | less 

awk '{OFS="\t"; print $2,$3,$4,$5,$7,$9,$11,$12}' SNP150_common.txt > SNP150_table.txt
awk 'NR>1 {OFS="\t"; print $2,$3-1,$4+1,$5}' SNP150_common.txt | awk 'NR < 1489644 {OFS="\t"; if ( $2 < 0 ) print $1,"0",$3,$4; else print $0}' > SNP150_common.bed


sort -k 1,1 -k2,2n eQTL_peakIM_check_vars2.txt > eQTL_peakIM_check_vars2.bed
sort -k 1,1 -k2,2n SNP150.bed > SNP150_sort.bed
grep "rs61960554" SNP150_sort.bed 

awk '{print $1}' SNP150_sort.bed | uniq
module load bedtools2
bedtools intersect -a  eQTL_peakIM_check_vars2.bed -b SNP150_common.bed -wo | less 
bedtools intersect -a  eQTL_sumIM_check_vars2.txt -b SNP150.bed -wo | less 


3. Could try to get an EBV estimate using the bisulphite reads.. ?
done.

How does the correlation match up?



5. correlating the ATAC_signal at summits with annotated gene expression?
make directory
```bash
cd /gs/gsfs0/users/anjohnst/17_member/QTLs_2017/Integration
mkdir ATAC_Gene_Cor
cd ATAC_Gene_Cor
module load R/3.4.0
```

```

```


