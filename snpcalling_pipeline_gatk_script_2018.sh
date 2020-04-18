#!/bin/sh

###################################################################
# pipeline for running snp-calling pipeline on re-sequencing data
#
# usage:  snpcalling_pipeline_gatk_script_2018.sh $basedir $num_of_processors
#
###################################################################

cd $1
CPU=$2

hisat2-build ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00.fa ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00

for file in `dir -d *_R1_001.fastq ` ; do
    file2=`echo "$file" | sed 's/_R1_001.fastq/_R2_001.fastq/'`
    sortedbamfile=`echo "$file" | sed 's/_R1_001.fastq/.sort.bam/'`
#    samfile=`echo "$file" | sed 's/.fastq/.sam/'`
    hisat2 -p $CPU -x ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00 -1 $file -2 $file2 | samtools view -bS -@ $CPU | samtools sort -@ $CPU - > $sortedbamfile
done


ls *.sort.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/programs/picard-tools-1.119/MarkDuplicates.jar INPUT={} OUTPUT={.}.md.bam METRICS_FILE={.}.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

ls *.md.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/programs/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT={} OUTPUT={}.rg.bam SORT_ORDER=coordinate RGID={} RGLB=1 RGPL=illumina RGPU=run RGSM={} RGCN=tom360 RGDS={}

ls *.rg.bam | parallel -j $CPU samtools index {}

ls *.rg.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/programs/gatk_3.4/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00.fa -U ALLOW_N_CIGAR_READS -nt 1 -I {}  -o {.}.intervals

ls *rg.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/programs/gatk_3.4/GenomeAnalysisTK.jar -T IndelRealigner -R ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00.fa -U ALLOW_N_CIGAR_READS -I {}  -o {.}.ralgn.bam -targetIntervals {.}.intervals

# -Djava.io.tmpdir=`pwd`~/tmp add this to the next command?:
ls *.ralgn.bam |parallel -j $CPU java -jar ~/programs/gatk_3.4/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00.fa -U ALLOW_N_CIGAR_READS -I {} -o {.}.g.vcf --emitRefConfidence GVCF

# ls *.ralgn.bam |parallel -j $CPU java -jar ~/programs/gatk_3.4/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00.fa  -I {} -o {.}.vcf -U ALLOW_N_CIGAR_READS

# GenotypeGVCFs
ls *vcf > lov.list
java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R ~/sequence_data/reference/S_lycopersicum_chromosomes.3.00 -V lov.list -nt $CPU -o combined.vcf