grep -E "MI0001064|MI0001065|MI0003726|MI0001067|MI0004988|MI0001068" /home/greally-lab/indexes/hg38/miRNA/ebv.gtf | less 
ebv-mir-BHRF1-1
ebv-miR-BHRF1-2-3p
ebv-miR-BART4-5p
ebv-miR-BART1-5p
ebv-miR-BART2-5p
ebv-miR-BART15 (overlap with premiRNA)

"while in the autologous LCL, of 44 miRNA, only a part of them were amplified (EBVmiR-BHRF1-1, miR-BHRF1-2 and miR-BHRF1-3 of the BHRF cluster, ebv-BART1-3p, BART1-5p, BART2-5p, BART3, BART3*, BART4 and BART15 of the BART cluster) (Table S1). These results are consistent with the deletion of 12 Kb in the BamHI-A region carried from the B95.8 strain [70]."

grep -E "MIMAT0000995|MIMAT0000997|MIMAT0003412|MIMAT0000999|MIMAT0001000" /home/greally-lab/indexes/hg38/miRNA/ebv.gtf | less 


chrEBV  .       miRNA_primary_transcript        41471   41536   .       +       .       ID "MI0001064"; Alias "MI0001064"; Name "ebv-mir-BHRF1-1";
chrEBV  .       miRNA   41474   41495   .       +       .       ID "MIMAT0000995"; Alias "MIMAT0000995"; Name "ebv-miR-BHRF1-1"; Derives_from "MI0001064"
chrEBV  .       miRNA_primary_transcript        42848   42912   .       +       .       ID "MI0001065"; Alias "MI0001065"; Name "ebv-mir-BHRF1-2";
chrEBV  .       miRNA   42853   42874   .       +       .       ID "MIMAT0000996"; Alias "MIMAT0000996"; Name "ebv-miR-BHRF1-2-5p"; Derives_from "MI0001065"
chrEBV  .       miRNA   42888   42909   .       +       .       ID "MIMAT0000997"; Alias "MIMAT0000997"; Name "ebv-miR-BHRF1-2-3p"; Derives_from "MI0001065"
chrEBV  .       miRNA_primary_transcript        139220  139295  .       +       .       ID "MI0003726"; Alias "MI0003726"; Name "ebv-mir-BART4";
chrEBV  .       miRNA   139228  139249  .       +       .       ID "MIMAT0003412"; Alias "MIMAT0003412"; Name "ebv-miR-BART4-5p"; Derives_from "MI0003726"
chrEBV  .       miRNA   139266  139288  .       +       .       ID "MIMAT0009204"; Alias "MIMAT0009204"; Name "ebv-miR-BART4-3p"; Derives_from "MI0003726"
chrEBV  .       miRNA_primary_transcript        139346  139415  .       +       .       ID "MI0001067"; Alias "MI0001067"; Name "ebv-mir-BART1";
chrEBV  .       miRNA   139351  139374  .       +       .       ID "MIMAT0000999"; Alias "MIMAT0000999"; Name "ebv-miR-BART1-5p"; Derives_from "MI0001067"
chrEBV  .       miRNA   139387  139408  .       +       .       ID "MIMAT0003390"; Alias "MIMAT0003390"; Name "ebv-miR-BART1-3p"; Derives_from "MI0001067"
chrEBV  .       miRNA_primary_transcript        139507  139584  .       +       .       ID "MI0004988"; Alias "MI0004988"; Name "ebv-mir-BART15";
chrEBV  .       miRNA   139553  139574  .       +       .       ID "MIMAT0003713"; Alias "MIMAT0003713"; Name "ebv-miR-BART15"; Derives_from "MI0004988"
chrEBV  .       miRNA_primary_transcript        152745  152806  .       +       .       ID "MI0001068"; Alias "MI0001068"; Name "ebv-mir-BART2";
chrEBV  .       miRNA   152747  152768  .       +       .       ID "MIMAT0001000"; Alias "MIMAT0001000"; Name "ebv-miR-BART2-5p"; Derives_from "MI0001068"
chrEBV  .       miRNA   152783  152806  .       +       .       ID "MIMAT0004744"; Alias "MIMAT0004744"; Name "ebv-miR-BART2-3p"; Derives_from "MI0001068"

chrEBV  .       miRNA   41474   41495   .       +       .       ID "MIMAT0000995"; Alias "MIMAT0000995"; Name "ebv-miR-BHRF1-1"; Derives_from "MI0001064"
chrEBV  .       miRNA   42888   42909   .       +       .       ID "MIMAT0000997"; Alias "MIMAT0000997"; Name "ebv-miR-BHRF1-2-3p"; Derives_from "MI0001065"
chrEBV  .       miRNA   139228  139249  .       +       .       ID "MIMAT0003412"; Alias "MIMAT0003412"; Name "ebv-miR-BART4-5p"; Derives_from "MI0003726"
chrEBV  .       miRNA   139351  139374  .       +       .       ID "MIMAT0000999"; Alias "MIMAT0000999"; Name "ebv-miR-BART1-5p"; Derives_from "MI0001067"
chrEBV  .       miRNA   152747  152768  .       +       .       ID "MIMAT0001000"; Alias "MIMAT0001000"; Name "ebv-miR-BART2-5p"; Derives_from "MI0001068"


What about quantifying the alignments the the EBV genome?

for f1 in ../*Aligned.sortedByCoord.out.bam;
do
SAMPLE="$(echo ${f1} | cut -d '/' -f2 | cut -d '.' -f1 | cut -c1-4)"
echo $SAMPLE
qsub -S /bin/bash -N ${SAMPLE}_featureCounts -j y -cwd -q highmem.q -pe smp-highmem 2 -l h_vmem=5.6G << EOF
module load subread/1.5.0-p1/gcc.4.4.7 
featureCounts -T 2 -F GTF -t miRNA -g ID -a /home/greally-lab/indexes/hg38/miRNA/ebv.gtf -o Count_miRNA_EBV/${SAMPLE}_EBV.txt ${f1}
EOF
done

grep -f EBV_express_names.txt /home/greally-lab/indexes/hg38/miRNA/ebv.gtf | less 

chrEBV  .       miRNA   41474   41495   .       +       .       ID "MIMAT0000995"; Alias "MIMAT0000995"; Name "ebv-miR-BHRF1-1"; Derives_from "MI0001064"
chrEBV  .       miRNA   42853   42874   .       +       .       ID "MIMAT0000996"; Alias "MIMAT0000996"; Name "ebv-miR-BHRF1-2-5p"; Derives_from "MI0001065"
chrEBV  .       miRNA   42888   42909   .       +       .       ID "MIMAT0000997"; Alias "MIMAT0000997"; Name "ebv-miR-BHRF1-2-3p"; Derives_from "MI0001065"
chrEBV  .       miRNA   42968   42989   .       +       .       ID "MIMAT0000998"; Alias "MIMAT0000998"; Name "ebv-miR-BHRF1-3"; Derives_from "MI0001066"
chrEBV  .       miRNA   139087  139107  .       +       .       ID "MIMAT0003410"; Alias "MIMAT0003410"; Name "ebv-miR-BART3-5p"; Derives_from "MI0003725"
chrEBV  .       miRNA   139124  139145  .       +       .       ID "MIMAT0003411"; Alias "MIMAT0003411"; Name "ebv-miR-BART3-3p"; Derives_from "MI0003725"
chrEBV  .       miRNA   139228  139249  .       +       .       ID "MIMAT0003412"; Alias "MIMAT0003412"; Name "ebv-miR-BART4-5p"; Derives_from "MI0003726"
chrEBV  .       miRNA   139266  139288  .       +       .       ID "MIMAT0009204"; Alias "MIMAT0009204"; Name "ebv-miR-BART4-3p"; Derives_from "MI0003726"
chrEBV  .       miRNA   139351  139374  .       +       .       ID "MIMAT0000999"; Alias "MIMAT0000999"; Name "ebv-miR-BART1-5p"; Derives_from "MI0001067"
chrEBV  .       miRNA   139387  139408  .       +       .       ID "MIMAT0003390"; Alias "MIMAT0003390"; Name "ebv-miR-BART1-3p"; Derives_from "MI0001067"
chrEBV  .       miRNA   139553  139574  .       +       .       ID "MIMAT0003713"; Alias "MIMAT0003713"; Name "ebv-miR-BART15"; Derives_from "MI0004988"
chrEBV  .       miRNA   152747  152768  .       +       .       ID "MIMAT0001000"; Alias "MIMAT0001000"; Name "ebv-miR-BART2-5p"; Derives_from "MI0001068"
chrEBV  .       miRNA   152783  152806  .       +       .       ID "MIMAT0004744"; Alias "MIMAT0004744"; Name "ebv-miR-BART2-3p"; Derives_from "MI0001068"


ebv-miR-BHRF1-1
ebv-miR-BHRF1-2-5p
ebv-miR-BHRF1-2-3p
ebv-miR-BHRF1-3

ebv-miR-BART3-5p
ebv-miR-BART3-3p

(EBVmiR-BHRF1-1, miR-BHRF1-2 and miR-BHRF1-3 of the BHRF cluster, ebv-BART1-3p, BART1-5p, BART2-5p, BART3, BART3*, BART4 and BART15 of the BART cluster) 








