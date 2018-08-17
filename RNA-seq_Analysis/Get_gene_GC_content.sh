module load bedtools2
bedtools nuc -fi /home/greally-lab/indexes/hg38_GenCode_EBV_decoy/Gencode_EBV_decoy_noComments.fasta -bed gencode.v24.primary_assembly.annotation.EBV.GENES.sorted.bed > Gene_GCContent.txt

awk '{print $4"\t"$8"\t"$15}' Gene_GCContent.txt > Gene_GCContent_DF.txt

