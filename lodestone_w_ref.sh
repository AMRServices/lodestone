#!/bin/bash

#########
### usage 
#########
usage()
{
cat << EOF
usage: sh $0 options

This script is the SP3 pipeline implementation
as shared by Stephen Bush

Author: Hardik Parikh
Last Modified: 01/21/2019
Version: 0.2

OPTIONS:
        -h      Show this message
        -c      CAVID
        -1      FULL/PATH/TO/ForwardRawReads (*_R1.fq.gz)
        -2      FULL/PATH/TO/ReverseRawReads (*_R2.fq.gz)
        -r      FULL/PATH/TO/ref.fa (Reference genome FASTA)
        -o      FULL/PATH/TO/OutDir
EOF
}

if [ $# -eq 0 ]; then
        usage
        exit
fi


##################
### Accept options
##################
while getopts “hc:1:2:r:o:?” OPTION
do
        case $OPTION in
                h)      usage
                        exit 1;;
                c)      CAVID=$OPTARG;;
                1)      IN_R1FILE=$OPTARG;;
                2)      IN_R2FILE=$OPTARG;;
                r)      REFFA=$OPTARG;;
                o)      OUTDIR=$OPTARG;;
                ?)      usage
                        exit 1;;

        esac
done


#####################
### Default variables
#####################
THREADS=20
LOGFILE=${OUTDIR}/log.${CAVID}.txt
R1FILE=${OUTDIR}/${CAVID}_R1.fq.gz
R2FILE=${OUTDIR}/${CAVID}_R2.fq.gz


##########
### Set up
##########
mkdir -p ${OUTDIR}
cd ${OUTDIR}
touch ${LOGFILE}

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. beginning analysis of ${CAVID}\n" >> ${LOGFILE}
cp ${IN_R1FILE} ${R1FILE}
cp ${IN_R2FILE} ${R2FILE}
cp ${REFFA} ref.fa

################
### Dependencies
################
# fastqc, trim_galore, trimmomatic
# bwa, picard, samtools, bcftools, vcftools
# spades, quast


###############
### Sequence QC
###############

## FastQC raw-data
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running initial FastQC on ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/1_initial_fastqc_report
fastqc ${R1FILE} ${R2FILE} &>> ${LOGFILE}
mv ${OUTDIR}/${CAVID}_R1_fastqc.html ${OUTDIR}/logs/1_initial_fastqc_report
mv ${OUTDIR}/${CAVID}_R2_fastqc.html ${OUTDIR}/logs/1_initial_fastqc_report
mv ${OUTDIR}/${CAVID}_R1_fastqc.zip ${OUTDIR}/logs/1_initial_fastqc_report
mv ${OUTDIR}/${CAVID}_R2_fastqc.zip ${OUTDIR}/logs/1_initial_fastqc_report

## TrimGalore adapter removal and trimming
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running Trim Galore on ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/2_trimgalore_report
trim_galore --stringency 5 --trim-n --max_n 5 --length 35 --paired -o ${OUTDIR} ${R1FILE} ${R2FILE} &>> ${LOGFILE}
mv ${R1FILE}_trimming_report.txt ${OUTDIR}/logs/2_trimgalore_report
mv ${R2FILE}_trimming_report.txt ${OUTDIR}/logs/2_trimgalore_report
mv ${OUTDIR}/${CAVID}_R1_val_1.fq.gz ${R1FILE}
mv ${OUTDIR}/${CAVID}_R2_val_2.fq.gz ${R2FILE}

## Trimmomatic read trimming
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running Trimmomatic on Trim Galore validated reads for ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/3_trimmomatic_report
java -jar $TRIMMOMATICPATH/trimmomatic-0.36.jar PE -phred33 ${R1FILE} ${R2FILE} ${OUTDIR}/${CAVID}.R1_paired.fq.gz ${OUTDIR}/${CAVID}.R1_unpaired.fq.gz ${OUTDIR}/${CAVID}.R2_paired.fq.gz ${OUTDIR}/${CAVID}.R2_unpaired.fq.gz TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70 > ${OUTDIR}/logs/3_trimmomatic_report/${CAVID}.trimmomatic_log.txt 2>&1
mv ${OUTDIR}/${CAVID}.R1_paired.fq.gz ${R1FILE}
mv ${OUTDIR}/${CAVID}.R2_paired.fq.gz ${R2FILE}
cat ${OUTDIR}/${CAVID}.R1_unpaired.fq.gz ${OUTDIR}/${CAVID}.R2_unpaired.fq.gz > ${OUTDIR}/${CAVID}.se.fq.gz
rm -rf ${OUTDIR}/${CAVID}.R1_unpaired.fq.gz ${OUTDIR}/${CAVID}.R2_unpaired.fq.gz

## FastQC hq reads
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running FastQC on high-quality reads for ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/4_final_fastqc_report
fastqc ${R1FILE} ${R2FILE} ${OUTDIR}/${CAVID}.se.fq.gz &>> ${LOGFILE}
mv ${OUTDIR}/${CAVID}_R1_fastqc.html ${OUTDIR}/logs/4_final_fastqc_report
mv ${OUTDIR}/${CAVID}_R2_fastqc.html ${OUTDIR}/logs/4_final_fastqc_report
mv ${OUTDIR}/${CAVID}.se_fastqc.html ${OUTDIR}/logs/4_final_fastqc_report
mv ${OUTDIR}/${CAVID}_R1_fastqc.zip ${OUTDIR}/logs/4_final_fastqc_report
mv ${OUTDIR}/${CAVID}_R2_fastqc.zip ${OUTDIR}/logs/4_final_fastqc_report
mv ${OUTDIR}/${CAVID}.se_fastqc.zip ${OUTDIR}/logs/4_final_fastqc_report


######################
### Align Reads to Ref
######################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. aligning reads to reference genome for ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/8_bam_statistics

## Create BWA Index for reference
bwa index ${OUTDIR}/ref.fa 2>> ${LOGFILE}

## Align PE reads
bwa mem -R "@RG\tID:group_${CAVID}\tSM:sample_${CAVID}\tPL:Illumina\tLIB:lib_${CAVID}\tPU:unit_${CAVID}" -t ${THREADS} -M ${OUTDIR}/ref.fa ${R1FILE} ${R2FILE} 2>> ${LOGFILE} | sambamba view -S -f bam -t ${THREADS} -o ${OUTDIR}/${CAVID}.unsorted.pe.bam /dev/stdin
java -jar $PICARD CleanSam INPUT=${OUTDIR}/${CAVID}.unsorted.pe.bam OUTPUT=${OUTDIR}/${CAVID}.unsorted.cleaned.pe.bam TMP_DIR=${OUTDIR} &>> ${LOGFILE}
java -jar $PICARD FixMateInformation INPUT=${OUTDIR}/${CAVID}.unsorted.cleaned.pe.bam OUTPUT=${OUTDIR}/${CAVID}.unsorted.cleaned.fixedmate.pe.bam TMP_DIR=${OUTDIR} &>> ${LOGFILE}
java -jar $PICARD SortSam INPUT=${OUTDIR}/${CAVID}.unsorted.cleaned.fixedmate.pe.bam OUTPUT=${OUTDIR}/${CAVID}.sorted.pe.bam SORT_ORDER=coordinate TMP_DIR=${OUTDIR} &>> ${LOGFILE}
java -jar $PICARD MarkDuplicates INPUT=${OUTDIR}/${CAVID}.sorted.pe.bam OUTPUT=${OUTDIR}/${CAVID}.pe.bam METRICS_FILE=${OUTDIR}/${CAVID}.pe.metrics ASSUME_SORTED=true TMP_DIR=${OUTDIR} &>> ${LOGFILE}
samtools view -b -F 2048 -F 256 ${OUTDIR}/${CAVID}.pe.bam > ${OUTDIR}/${CAVID}.pe.rm_supplementary.bam 2>> ${LOGFILE}
mv ${OUTDIR}/${CAVID}.pe.rm_supplementary.bam ${OUTDIR}/${CAVID}.pe.bam
java -jar $PICARD BuildBamIndex INPUT=${OUTDIR}/${CAVID}.pe.bam &>> ${LOGFILE}

## Align SE reads
bwa mem -R "@RG\tID:group_${CAVID}\tSM:sample_${CAVID}\tPL:Illumina\tLIB:lib_${CAVID}\tPU:unit_${CAVID}" -t ${THREADS} -M ${OUTDIR}/ref.fa ${OUTDIR}/${CAVID}.se.fq.gz 2>> ${LOGFILE} | sambamba view -S -f bam -t ${THREADS} -o ${OUTDIR}/${CAVID}.unsorted.se.bam /dev/stdin
java -jar $PICARD CleanSam INPUT=${OUTDIR}/${CAVID}.unsorted.se.bam OUTPUT=${OUTDIR}/${CAVID}.unsorted.cleaned.se.bam TMP_DIR=${OUTDIR} &>> ${LOGFILE}
java -jar $PICARD SortSam INPUT=${OUTDIR}/${CAVID}.unsorted.cleaned.se.bam OUTPUT=${OUTDIR}/${CAVID}.sorted.se.bam SORT_ORDER=coordinate TMP_DIR=${OUTDIR} &>> ${LOGFILE}
java -jar $PICARD MarkDuplicates INPUT=${OUTDIR}/${CAVID}.sorted.se.bam OUTPUT=${OUTDIR}/${CAVID}.se.bam METRICS_FILE=${OUTDIR}/${CAVID}.se.metrics ASSUME_SORTED=true TMP_DIR=${OUTDIR} &>> ${LOGFILE}
samtools view -b -F 2048 -F 256 ${OUTDIR}/${CAVID}.se.bam > ${OUTDIR}/${CAVID}.se.rm_supplementary.bam 2>> ${LOGFILE}
mv ${OUTDIR}/${CAVID}.se.rm_supplementary.bam ${OUTDIR}/${CAVID}.se.bam
java -jar $PICARD BuildBamIndex INPUT=${OUTDIR}/${CAVID}.se.bam &>> ${LOGFILE}


## remove intermediate files
rm ${OUTDIR}/${CAVID}.unsorted.pe.bam ${OUTDIR}/${CAVID}.unsorted.cleaned.pe.bam ${OUTDIR}/${CAVID}.unsorted.cleaned.fixedmate.pe.bam ${CAVID}.sorted.pe.bam 
rm ${OUTDIR}/${CAVID}.unsorted.se.bam ${OUTDIR}/${CAVID}.unsorted.cleaned.se.bam ${CAVID}.sorted.se.bam 

## Merge BAMS
java -jar $PICARD MergeSamFiles SORT_ORDER='coordinate' INPUT=${OUTDIR}/${CAVID}.pe.bam INPUT=${OUTDIR}/${CAVID}.se.bam OUTPUT=${OUTDIR}/${CAVID}.bam TMP_DIR=${OUTDIR} CREATE_INDEX=TRUE &>> ${LOGFILE}

## Get BAM Stats
samtools stats ${OUTDIR}/${CAVID}.bam > ${OUTDIR}/logs/8_bam_statistics/${CAVID}.samtools_stats_output.txt
samtools flagstat ${OUTDIR}/${CAVID}.bam > ${OUTDIR}/logs/8_bam_statistics/${CAVID}.samtools_flagstat_output.txt


#################
### Call Variants
#################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. calling variants against reference genome for ${CAVID}...\n" >> ${LOGFILE}

bcftools mpileup -Ou -f ${OUTDIR}/ref.fa ${OUTDIR}/${CAVID}.bam 2>> ${LOGFILE} | bcftools call --threads ${THREADS} --ploidy 1 -mv -Ov -o ${OUTDIR}/${CAVID}.vcf
vcfallelicprimitives ${OUTDIR}/${CAVID}.vcf > ${OUTDIR}/${CAVID}.regularised.vcf 2>> ${LOGFILE}
mv ${OUTDIR}/${CAVID}.regularised.vcf ${OUTDIR}/${CAVID}.vcf
vcf-sort ${OUTDIR}/${CAVID}.vcf > ${OUTDIR}/${CAVID}.sorted.vcf 2>> ${LOGFILE}
mv ${OUTDIR}/${CAVID}.sorted.vcf ${OUTDIR}/${CAVID}.vcf
/apps/software/standard/core/htslib/1.9/bin/bgzip ${OUTDIR}/${CAVID}.vcf
/apps/software/standard/core/htslib/1.9/bin/tabix -p vcf ${OUTDIR}/${CAVID}.vcf.gz


###########################
### SPAdes de-novo assembly
###########################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. SPAdes de-novo assembly for ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/9_assembly_statistics
spades.py --cov-cutoff auto --careful -1 ${R1FILE} -2 ${R2FILE} -s ${OUTDIR}/${CAVID}.se.fq.gz -t 20 -o ${OUTDIR}/spades &>> ${LOGFILE}
mv ${OUTDIR}/spades/contigs.fasta ${OUTDIR}/${CAVID}.fa
### check this ###
quast.py -t ${THREADS} --pe1 ${R1FILE} --pe2 ${R2FILE} --single ${OUTDIR}/${CAVID}.se.fq.gz -o ${OUTDIR}/logs/9_assembly_statistics/quast_out ${OUTDIR}/${CAVID}.fa &>> ${LOGFILE}


############
### Clean up
############
## remove intermediate files
rm ${OUTDIR}/ftpdirpaths.txt
rm ${OUTDIR}/${CAVID}*metrics
rm ${OUTDIR}/${CAVID}.pe.bam ${OUTDIR}/${CAVID}.pe.bai ${OUTDIR}/${CAVID}.se.bam ${OUTDIR}/${CAVID}.se.bai
md5sum ${OUTDIR}/${CAVID}.fa > ${OUTDIR}/${CAVID}.fa.md5sum
md5sum ${OUTDIR}/${CAVID}.vcf.gz > ${OUTDIR}/${CAVID}.vcf.gz.md5sum
md5sum ${OUTDIR}/${CAVID}.vcf.gz.tbi > ${OUTDIR}/${CAVID}.vcf.gz.tbi.md5sum


############
### DONE ###
############
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. Finished processing ${CAVID}!\n" >> ${LOGFILE}
exit 0
