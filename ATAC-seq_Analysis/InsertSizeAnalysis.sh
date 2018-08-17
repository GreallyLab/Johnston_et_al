for s in Bams/*_UMap_noMT_mkdup.bam
do
dest1="${s##*/}"
echo $dest1
SampleName="$(echo ${dest1} | cut -d '_' -f1)"
echo $SampleName
DIR='InsertSizeAnalysis/'
DIR_ANNO='/home/ajohnsto/17_member/ATACseq/Replicates/Annotations/'
	qsub -S /bin/bash -N ${SampleName}_Insert_Analysis -cwd -l h_vmem=30G -j y << EOF
	module load samtools 
	module load picard-tools/1.92/java.1.8.0_20
	module load bedtools2
	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_InsertSizeMetrics.pdf INPUT=${s} OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_InsertSizeMetrics.txt ASSUME_SORTED=false DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false

 	intersectBed -a ${s} -b $DIR_ANNO/hg19-ActivePromoter_hg38lift.bed > $DIR${SampleName}_AP.bam
 	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_AP.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_AP.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_AP.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
 	
 	intersectBed -a ${s} -b $DIR_ANNO/hg19-WeakPromoter_hg38lift.bed > $DIR${SampleName}_WP.bam
	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_WP.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_WP.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_WP.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
 	
 	intersectBed -a ${s} -b $DIR_ANNO/hg19-StrongEnhancer_hg38lift.bed > $DIR${SampleName}_SE.bam
	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_SE.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_SE.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_SE.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
 	
 	intersectBed -a ${s} -b $DIR_ANNO/hg19-WeakEnhancer_hg38lift.bed > $DIR${SampleName}_WE.bam
 	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_WE.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_WE.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_WE.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
 	
 	intersectBed -a ${s} -b $DIR_ANNO/hg19-TxnAll_hg38lift.bed > $DIR${SampleName}_Txn.bam
	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_Txn.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_Txn.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_Txn.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
 	
 	intersectBed -a ${s} -b $DIR_ANNO/hg19-Repressed_hg38lift.bed > $DIR${SampleName}_Repressed.bam
	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_Repressed.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_Repressed.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_Repressed.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
 	
 	intersectBed -a ${s} -b $DIR_ANNO/hg19-Heterochrom_hg38lift.bed > $DIR${SampleName}_Heterochrom.bam
	java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_Heterochrom.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_Heterochrom.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_Heterochrom.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false

intersectBed -a ${s} -b $DIR_ANNO/hg19-refSeq-promTSS_hg38lift.bed > $DIR${SampleName}_promTSS.bam
java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_promTSS.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_promTSS.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_promTSS.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false	

module load R
	
Rscript /home/ajohnsto/17_member/ATACseq/Rep1/Fastq/comb_reads/InsertSizeAnalysis/InsertSizeAnalysis.R $DIR${SampleName}
EOF
done

less $DIR_ANNO/hg19-refSeq-promTSS_hg38lift.bed

	66533955	66751139

	Rscript /home/ajohnsto/17_member/ATACseq/Rep1/Fastq/comb_reads/InsertSizeAnalysis/InsertSizeAnalysis.R $DIR${SampleName}


cannot open file 'InsertSizeAnalysis/gm87new_UMap_rmChrM_rmdup_InsertSizeMetrics.txt': No such file or directory

gm87new_UMap_noMT_mkdup.bam

intersectBed -a Bams/gm87new_UMap_noMT_mkdup.bam -b $DIR_ANNO/hg19-refSeq-promTSS_hg38lift.bed > $DIR${SampleName}_promTSS.bam

for s in Bams/gm87new_UMap_noMT_mkdup.bam
do
dest1="${s##*/}"
echo $dest1
SampleName="$(echo ${dest1} | cut -d '_' -f1)"
echo $SampleName
DIR='InsertSizeAnalysis/'
DIR_ANNO='/home/ajohnsto/17_member/ATACseq/Replicates/Annotations/'
	qsub -S /bin/bash -N ${SampleName}_Insert_Analysis -cwd -l h_vmem=30G -j y << EOF
	module load samtools 
	module load picard-tools/1.92/java.1.8.0_20
	module load bedtools2
intersectBed -a ${s} -b $DIR_ANNO/hg19-refSeq-promTSS_hg38lift.bed > $DIR${SampleName}_promTSS.bam
java -jar $(which CollectInsertSizeMetrics.jar) HISTOGRAM_FILE=$DIR${SampleName}_UMap_mkdup_noMT_promTSS.InsertSizeMetrics.pdf INPUT=$DIR${SampleName}_promTSS.bam OUTPUT=$DIR${SampleName}_UMap_mkdup_noMT_promTSS.InsertSizeMetrics.txt ASSUME_SORTED=true DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=ALL_READS STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false	
module load R
Rscript /home/ajohnsto/17_member/ATACseq/Rep1/Fastq/comb_reads/InsertSizeAnalysis/InsertSizeAnalysis.R $DIR${SampleName}
EOF
done
