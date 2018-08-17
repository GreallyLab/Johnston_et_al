# Redoing the peak analysis using +/- 250bp around peaks summits

qsub -S /bin/bash -N Rep1_call_summit -cwd -l h_vmem=30G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
macs2 callpeak --nomodel --nolambda -g 3e9 --keep-dup 'all' --call-summits -t Rep1_shifted_sort.bed -n ../Peaks_Rep1/All_Rep1_sort
EOF

qsub -S /bin/bash -N Rep2_call_summit -cwd -l h_vmem=30G -j y << EOF
module load MACS2/2.1.0/python.2.7.8
macs2 callpeak --nomodel --nolambda -g 3e9 --keep-dup 'all' --call-summits -t Rep2_shifted.bed -n ../Peaks_Rep2/All_Rep2_sort
EOF


# ENCODE mappability consensus blacklist downloaded.
# wgEncodeDacMapabilityConsensusExcludable.bed.gz from 
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/
# to /home/greally-lab/indexes

# check out http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/
# downloaded to to /home/greally-lab/indexes
wc -l Kundeje_hg38_blackList.bed
# 38 Kundeje_hg38_blackList.bed
module load bedtools2 
bedtools intersect -a Kundeje_hg38_blackList.bed -b consensusBlacklist_hg38.bed  -u | wc -l # 3

awk '{print $1}' Kundeje_hg38_blackList.bed | uniq -c # 10, 16, 1, 20, 21, 2, 3, 4, 5
awk '{print $1}' consensusBlacklist_hg38.bed | uniq -c

bedtools intersect -a Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -b /home/greally-lab/indexes/Kundeje_hg38_blackList.bed -u | wc -l
# 31, that's odd

bedtools intersect -a Summit_idrPeak_chr_250_noOver_sort_merge_noOver_Canon_FilMappMitoBlack_dupfix_sort.bed -b /home/greally-lab/indexes/Kundeje_hg38_blackList.bed -wo | less

wc -l wgEncodeHg19ConsensusSignalArtifactRegions.bed
411 wgEncodeHg19ConsensusSignalArtifactRegions.bed

bedtools intersect -a wgEncodeHg19ConsensusSignalArtifactRegions.bed -b consensusBlacklist.bed  -u | wc -l

# it's in hg19 so I need to liftover to hg38

~/Programs/LiftOverProg/liftOver consensusBlacklist.bed ~/Programs/LiftOverProg/hg19ToHg38.over.chain "consensusBlacklist_hg38.bed" "consensusBlacklist_hg38_unmapped.bed"

wc -l consensusBlacklist.bed
# 411 consensusBlacklist.bed

wc -l consensusBlacklist_hg38.bed
# 401 consensusBlacklist_hg38.bed

wc -l consensusBlacklist_hg38_unmapped.bed
# 20 consensusBlacklist_hg38_unmapped.bed, either split or partially deleted in the new one. 

# need to add +/-250 to each peak summit 

awk '{print $1"\t"$2-250"\t"$3+250"\t"$4"\t"$5"\t".}' File


# turns out that the IDR analysis that I did before has summits in the tenth column
# in  ~/17_member/ATACseq/Total_peak_analysis

# Reps_peaks_idr_05_chr

# bringing up the RProj which is where I made the peaks.bed file 

# Made a bed file for summits : Rep_nonspecific_summit.bed
awk '{print $1}' Rep_nonspecific_summit.bed | sort | uniq -c

# let's remove the blacklisted regions as well as the no canon chr # see bed Rep_nonspecific_summit_Canon.bed

# then remove the mappability blacklist regions: 
/home/greally-lab/indexes/consensusBlacklist_hg38.bed

module load bedtools2 
bedtools intersect -a Rep_nonspecific_summit_Canon.bed -b /home/greally-lab/indexes/consensusBlacklist_hg38.bed -v > Rep_nonspecific_summit_Canon_FilMappBlack.bed

wc -l Rep_nonspecific_summit_Canon.bed 
# 101741 Rep_nonspecific_summit_Canon.bed

wc -l Rep_nonspecific_summit_Canon_FilMappBlack.bed
# 101530 Rep_nonspecific_summit_Canon_FilMappBlack.bed

# removed 211 summit regions

# Then remove the mito homologue blacklisted regions: Mito_peak_blacklist.bed


module load bedtools2 
bedtools intersect -a Rep_nonspecific_summit_Canon_FilMappBlack.bed -b Mito_peak_blacklist.bed -v > Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed

wc -l Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed
# 101527 Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed
# removed 3 peaks that were mitochondrial homologues 

# let's get the quanitifcations 

# Rep1 
for f1 in ~/17_member/ATACseq/Rep1/Fastq/comb_reads/Beds/*_shifted_chr_noNeg.bed 
#for f1 in ~/17_member/ATACseq/Rep1/Fastq/comb_reads/Beds/gm77_shifted_chr_noNeg.bed 
do
SampleName="$(echo ${f1} | cut -d '/' -f12 | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N Quant_summit_${SampleName} -cwd -l h_vmem=10G -j y << EOF
module load bedtools2
bedtools intersect -a ${f1} -b ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed -wo | awk '{print \$10}' | cut -d '_' -f2 | sort -n | uniq -c > ${SampleName}_rep1_readPeaks
EOF
done

for f1 in ~/17_member/ATACseq/Replicates/FastQ/Beds/*_shifted_chr_noNeg.bed 
#for f1 in ~/17_member/ATACseq/Replicates/FastQ/Beds/gm77_shifted_chr_noNeg.bed 
do
SampleName="$(echo ${f1} | cut -d '/' -f11 | cut -d '_' -f1)"
echo $SampleName
qsub -S /bin/bash -N Quant_summit_${SampleName} -cwd -l h_vmem=10G -j y << EOF
module load bedtools2
bedtools intersect -a ${f1} -b ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed -wo | awk '{print \$10}' | cut -d '_' -f2 | sort -n | uniq -c > ${SampleName}_rep2_readPeaks
EOF
done


# For the CQN normalization, I need to have the GC content and length for each gene.
module load bedtools2
bedtools nuc -fi /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fasta -bed Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed > Rep_nonspecific_summit_Canon_FilMappMitoBlack_GCContent.txt

awk '{print $4"\t"$8"\t"$15}' Rep_nonspecific_summit_Canon_FilMappMitoBlack_GCContent.txt > Rep_nonspecific_summit_Canon_FilMappMitoBlack_GCContent_DF.txt

# ~/17_member/QTLs_2017 to get QTLs

cp ~/17_member/ATACseq/Total_peak_analysis/Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed .

module load bedtools2
bedtools intersect -a Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed -b haplotype_grch38_named_4intersect.txt -wo > ATAC_summit_Plat_haplo_interesect.txt

wc -l ATAC_summit_Plat_haplo_interesect.txt 
# 101,321 ATAC_summit_Plat_haplo_interesect.txt

bedtools intersect -a Rep_nonspecific_summit_Canon_FilMappMitoBlack.bed -b haplotype_grch38_named_4intersect.txt -wa | uniq -d
# no peaks overlap haplotypes


