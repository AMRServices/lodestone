#!/bin/bash
#$ -N SRR1582856
#$ -cwd
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-->$current_date_time. beginning analysis of SRR1582856
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "
-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. testing software availability before beginning..." >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### FastQC ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/fastqc/0.11.5/fastqc --version >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### cutadapt ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/cutadapt/1.16/bin/cutadapt --version >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### TrimGalore ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/trimgalore/0.4.5/trim_galore --version >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Trimmomatic ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /project/som_infc_dis_mathers/tools/Trimmomatic-0.36/trimmomatic-0.36.jar -version >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Kaiju ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaiju -h &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Mash ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/mash-Linux64-v2.1/bin/mash --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### seqtk ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/seqtk/bin/seqtk &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Entrez Direct ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/edirect/efetch --help &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### BWA ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/compiler/gcc/7.1.0/bwa/0.7.17/bin/bwa &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### SAMtools ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/samtools &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### BCFtools ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/bcftools --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Picard Tools ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar MarkDuplicates --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### VCFlib ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/vcflib/bin/vcfoverlay --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### VCFtools ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/compiler/gcc/7.1.0/vcftools/0.1.15/bin/vcf-sort --help &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### bgzip ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/bgzip --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### tabix ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/tabix --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### SPAdes ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
python /project/som_infc_dis_mathers/tools/SPAdes-3.11.1-Linux/bin/spades.py --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Quast ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
python /project/som_infc_dis_mathers/tools/quast-5.0.1/bin/quast.py --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Circos ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/circos-0.69-6/bin/circos --version &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### TETyper ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/genomics/runUvaBiopipe/uvabiopy/TEtyper.py -h &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "
-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. confirming database availability..." >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "
### Database employed by Kaiju ###
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
date_last_modified="`date -r /project/som_infc_dis_mathers/tools/kaiju/kaijudb/kaiju_db.fmi`"
echo "/project/som_infc_dis_mathers/tools/kaiju/kaijudb/kaiju_db.fmi. last modified $date_last_modified" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "
-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. beginning pipeline...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "the last command run - ${last_command} - finished with exit code $?."' EXIT
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. downloading fastqs for SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856
cd /scratch/hp7d/lodestone-testing/SRR1582856
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running initial FastQC on SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/fastqc/0.11.5/fastqc ./S1.fq.gz ./S2.fq.gz &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/1_initial_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1_fastqc.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/1_initial_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2_fastqc.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/1_initial_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1_fastqc.zip /scratch/hp7d/lodestone-testing/SRR1582856/logs/1_initial_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2_fastqc.zip /scratch/hp7d/lodestone-testing/SRR1582856/logs/1_initial_fastqc_report
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running cutadapt/TrimGalore on SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/trimgalore/0.4.5/trim_galore --path_to_cutadapt /apps/software/standard/core/cutadapt/1.16/bin/cutadapt --stringency 5 --trim-n --max_n 5 --paired --length 0 -o /scratch/hp7d/lodestone-testing/SRR1582856 ./S1.fq.gz ./S2.fq.gz &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/2_trimgalore_report
mv ./S1.fq.gz_trimming_report.txt /scratch/hp7d/lodestone-testing/SRR1582856/logs/2_trimgalore_report
mv ./S2.fq.gz_trimming_report.txt /scratch/hp7d/lodestone-testing/SRR1582856/logs/2_trimgalore_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1_val_1.fq.gz ./S1.fq.gz
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2_val_2.fq.gz ./S2.fq.gz
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/3_trimmomatic_report
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running Trimmomatic on the adapter-trimmed reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /project/som_infc_dis_mathers/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 20 ./S1.fq.gz ./S2.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.forward_paired.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.forward_unpaired.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.reverse_paired.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.reverse_unpaired.fq.gz TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100 >/scratch/hp7d/lodestone-testing/SRR1582856/logs/3_trimmomatic_report/SRR1582856.trimmomatic_log.txt 2>&1
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.forward_paired.fq.gz ./S1.fq.gz
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.reverse_paired.fq.gz ./S2.fq.gz
zcat /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.forward_unpaired.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.reverse_unpaired.fq.gz | gzip > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.forward_unpaired.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.reverse_unpaired.fq.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running FastQC on the cleaned and adapter-trimmed reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/fastqc/0.11.5/fastqc ./S1.fq.gz ./S2.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1_fastqc.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2_fastqc.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1_fastqc.zip /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2_fastqc.zip /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se_fastqc.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se_fastqc.zip /scratch/hp7d/lodestone-testing/SRR1582856/logs/4_final_fastqc_report
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running Kaiju on the cleaned and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaiju -z 20 -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -f /project/som_infc_dis_mathers/tools/kaiju/kaijudb/kaiju_db.fmi -i ./S1.fq.gz -j ./S2.fq.gz -a greedy -e 5 -E 0.05 -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.txt -r genus -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_genus_prediction.paired-end.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.txt -r species -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_species_prediction.paired-end.txt
awk '$1=="C" { print $2 }' /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.txt > /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.paired-end.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaiju2krona -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.txt -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report_for_krona.paired-end.txt &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.txt
/project/som_infc_dis_mathers/tools/KronaTools-2.7/bin/ktImportText -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.paired-end.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report_for_krona.paired-end.txt &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report_for_krona.paired-end.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running Kaiju on the cleaned and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaiju -z 20 -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -f /project/som_infc_dis_mathers/tools/kaiju/kaijudb/kaiju_db.fmi -i /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz -a greedy -e 5 -E 0.05 -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.txt -r genus -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_genus_prediction.single-end.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.txt -r species -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_species_prediction.single-end.txt
awk '$1=="C" { print $2 }' /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.txt > /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.single-end.txt
/project/som_infc_dis_mathers/tools/kaiju/bin/kaiju2krona -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.txt -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report_for_krona.single-end.txt &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.txt
/project/som_infc_dis_mathers/tools/KronaTools-2.7/bin/ktImportText -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report.single-end.html /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report_for_krona.single-end.txt &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_report_for_krona.single-end.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing reads unclassified by Kaiju from the cleaned and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/seqtk/bin/seqtk subseq ./S1.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.paired-end.txt | gzip > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1.excl_unclassifiable.fastq.gz
/project/som_infc_dis_mathers/tools/seqtk/bin/seqtk subseq ./S2.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.paired-end.txt | gzip > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2.excl_unclassifiable.fastq.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing reads unclassified by Kaiju from the cleaned and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/seqtk/bin/seqtk subseq /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.single-end.txt > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.excl_unclassifiable.fastq.gz
rm /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.paired-end.txt /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.classified_read_list.single-end.txt
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_1.excl_unclassifiable.fastq.gz ./S1.fq.gz
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856_2.excl_unclassifiable.fastq.gz ./S2.fq.gz
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.excl_unclassifiable.fastq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. downloading reference genome...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
top_species_hit_pe=$(awk -F '\t' 'NR==3 {print $3}' /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_species_prediction.paired-end.txt)
top_species_hit_se=$(awk -F '\t' 'NR==3 {print $3}' /scratch/hp7d/lodestone-testing/SRR1582856/logs/5_kaiju_report/SRR1582856.kaiju_species_prediction.single-end.txt)
if [ "$top_species_hit_pe" != "$top_species_hit_se" ]; then echo "ERROR: candidate species for SRR1582856 as predicted by top hit for PE reads ($top_species_hit_pe) != top hit for SE reads ($top_species_hit_se); unable to proceed with alignment
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt; fi
if [ "$top_species_hit_pe" != "$top_species_hit_se" ]; then exit 1; fi
GENOME=$(/project/som_infc_dis_mathers/tools/edirect/esearch -db genome -query "$top_species_hit_pe"[orgn] | /project/som_infc_dis_mathers/tools/edirect/efetch -format docsum | tee "/scratch/hp7d/lodestone-testing/SRR1582856/$top_species_hit_pe.genome.esearch.docsum")
ACC=`echo $GENOME | /project/som_infc_dis_mathers/tools/edirect/xtract -pattern DocumentSummary  -element Assembly_Accession`
RESULT=$(/project/som_infc_dis_mathers/tools/edirect/esearch -db assembly -query "$ACC" | /project/som_infc_dis_mathers/tools/edirect/efetch -format docsum | tee "/scratch/hp7d/lodestone-testing/SRR1582856/$top_species_hit_pe.assembly.esearch.docsum")
FTPP=`echo $RESULT | /project/som_infc_dis_mathers/tools/edirect/xtract -pattern DocumentSummary -element FtpPath_GenBank`
TAXID=`echo $RESULT | /project/som_infc_dis_mathers/tools/edirect/xtract -pattern DocumentSummary -element Taxid`
SPECIESTAXID=`echo $RESULT | /project/som_infc_dis_mathers/tools/edirect/xtract -pattern DocumentSummary -element SpeciesTaxid`
BASENAME=`basename $FTPP`
FTPPATHG=$FTPP/$BASENAME'_genomic.fna.gz'
FTPPATHGFF=$FTPP/$BASENAME'_genomic.gff.gz'
rm "/scratch/hp7d/lodestone-testing/SRR1582856/$top_species_hit_pe.genome.esearch.docsum" "/scratch/hp7d/lodestone-testing/SRR1582856/$top_species_hit_pe.assembly.esearch.docsum"
echo "-->species predicted: $top_species_hit_pe
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "-->NCBI taxonomy ID (species): $SPECIESTAXID
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded
echo "TOP KAIJU HIT (SPECIES): $top_species_hit_pe" >> /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.reference_genome_download_log.txt
echo "NCBI TAXON ID (FOR SPECIES): $SPECIESTAXID" >> /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.reference_genome_download_log.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. downloading the set of complete genomes corresponding to the best hit taxonomy ID...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
date_last_modified="`date -r assembly_summary.txt`"
echo "--> downloaded ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt. last modified $date_last_modified" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
awk -v speciestaxid=$SPECIESTAXID -F "\t" '$7==speciestaxid && $12=="Complete Genome" && $11=="latest"{print $20}' /scratch/hp7d/lodestone-testing/SRR1582856/assembly_summary.txt > /scratch/hp7d/lodestone-testing/SRR1582856/ftpdirpaths.txt
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' /scratch/hp7d/lodestone-testing/SRR1582856/ftpdirpaths.txt > /scratch/hp7d/lodestone-testing/SRR1582856/ftpfilepaths.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/ftpdirpaths.txt
number_of_complete_genomes="`wc -l < /scratch/hp7d/lodestone-testing/SRR1582856/ftpfilepaths.txt`"
echo "--> assembly_summary.txt parsed to identify $number_of_complete_genomes complete genomes for species taxon ID $SPECIESTAXID" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
echo "NUMBER OF COMPLETE GENOMES ASSIGNED THIS SPECIES TAXON ID, OBTAINED AFTER PARSING ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt (DOWNLOADED $date_last_modified): $number_of_complete_genomes" >> /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.reference_genome_download_log.txt
echo "GENOME URLS:" >> /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.reference_genome_download_log.txt
cat /scratch/hp7d/lodestone-testing/SRR1582856/ftpfilepaths.txt >> /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.reference_genome_download_log.txt
mkdir /scratch/hp7d/lodestone-testing/SRR1582856/temp_dir_genome_storage
wget -i /scratch/hp7d/lodestone-testing/SRR1582856/ftpfilepaths.txt -P /scratch/hp7d/lodestone-testing/SRR1582856/temp_dir_genome_storage &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/ftpfilepaths.txt
find /scratch/hp7d/lodestone-testing/SRR1582856/temp_dir_genome_storage -maxdepth 1 -type f > /scratch/hp7d/lodestone-testing/SRR1582856/fa_filepaths_for_mash_sketch.txt
cat ./S1.fq.gz ./S2.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.concatenated_reads.fq.gz
/project/som_infc_dis_mathers/tools/mash-Linux64-v2.1/bin/mash sketch -p 20 -l /scratch/hp7d/lodestone-testing/SRR1582856/fa_filepaths_for_mash_sketch.txt -o /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/mash-Linux64-v2.1/bin/mash sketch -p 20 -o /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fq /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.concatenated_reads.fq.gz &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.concatenated_reads.fq.gz
/project/som_infc_dis_mathers/tools/mash-Linux64-v2.1/bin/mash dist -p 20 /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fq.msh /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa.msh > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.mash-dist_output.tsv
closest_mash_hit=$(basename $(sort -g -k3 /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.mash-dist_output.tsv | head -1 | cut -f2))
closest_mash_hit_root=${closest_mash_hit%_genomic.fna.gz}
path_to_closest_mash_hit=$(grep $closest_mash_hit_root /scratch/hp7d/lodestone-testing/SRR1582856/assembly_summary.txt | head -1 | awk '{print $(NF)}')
FTPPATHG=$path_to_closest_mash_hit/$closest_mash_hit_root'_genomic.fna.gz'
FTPPATHGFF=$path_to_closest_mash_hit/$closest_mash_hit_root'_genomic.gff.gz'
wget $FTPPATHG -O /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.gz &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
gunzip /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.gz
wget $FTPPATHGFF -O /scratch/hp7d/lodestone-testing/SRR1582856/ref.gff.gz &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
gunzip /scratch/hp7d/lodestone-testing/SRR1582856/ref.gff.gz
echo "URL FOR CLOSEST MASH HIT: $FTPPATHG" >> /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.reference_genome_download_log.txt
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.mash-dist_output.tsv /scratch/hp7d/lodestone-testing/SRR1582856/logs/6_reference_genomes_downloaded/SRR1582856.mash-dist_output.tsv
rm /scratch/hp7d/lodestone-testing/SRR1582856/assembly_summary.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/fa_filepaths_for_mash_sketch.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fq.msh /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa.msh
rm -r /scratch/hp7d/lodestone-testing/SRR1582856/temp_dir_genome_storage
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. indexing reference genome...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/compiler/gcc/7.1.0/bwa/0.7.17/bin/bwa index /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. aligning the quality- and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/compiler/gcc/7.1.0/bwa/0.7.17/bin/bwa mem -R '@RG\tID:group_SRR1582856\tSM:sample_SRR1582856\tPL:Illumina\tLIB:lib_SRR1582856\tPU:unit_SRR1582856' -t 20 -M /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa ./S1.fq.gz ./S2.fq.gz 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt | /project/som_infc_dis_mathers/tools/sambamba-0.6.8/bin/sambamba view -S -f bam -t 20 -o /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam /dev/stdin
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. cleaning BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar CleanSam INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.cleaned.paired-end.bam TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.cleaned.paired-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. fixing mate information in the BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar FixMateInformation INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.fixmated.paired-end.bam TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.fixmated.paired-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. sorting BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar SortSam INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.paired-end.bam SORT_ORDER=coordinate TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. duplicate-marking BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar MarkDuplicates INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.paired-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bam METRICS_FILE=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.metrics ASSUME_SORTED=true TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing supplementary and non-primary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/samtools view -b -F 2048 -F 256 /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bam > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.rm_supplementary.bam 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.rm_supplementary.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bam
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar BuildBamIndex INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bam &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.paired-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.paired-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.metrics
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. aligning the quality- and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/compiler/gcc/7.1.0/bwa/0.7.17/bin/bwa mem -R '@RG\tID:group_SRR1582856\tSM:sample_SRR1582856\tPL:Illumina\tLIB:lib_SRR1582856\tPU:unit_SRR1582856' -t 20 -M /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt | /project/som_infc_dis_mathers/tools/sambamba-0.6.8/bin/sambamba view -S -f bam -t 20 -o /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam /dev/stdin
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. cleaning BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar CleanSam INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.cleaned.single-end.bam TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.cleaned.single-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. fixing mate information in the BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar FixMateInformation INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.fixmated.single-end.bam TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.fixmated.single-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. sorting BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar SortSam INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.single-end.bam SORT_ORDER=coordinate TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. duplicate-marking BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar MarkDuplicates INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.single-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bam METRICS_FILE=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.metrics ASSUME_SORTED=true TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing supplementary and non-primary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/samtools view -b -F 2048 -F 256 /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bam > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.rm_supplementary.bam 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.rm_supplementary.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bam
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar BuildBamIndex INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bam &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.unsorted.single-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.single-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.metrics
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. merging PE and SE BAMs for the aligned, classified, quality- and adapter-trimmed reads from SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar MergeSamFiles SORT_ORDER='coordinate' INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bam INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bam OUTPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bam TMP_DIR=/scratch/hp7d/lodestone-testing/SRR1582856 &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
java -jar /apps/software/standard/core/picard/2.18.5/picard.jar BuildBamIndex INPUT=/scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bam &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.paired-end.bai /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.single-end.bai
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/7_bam_statistics
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. obtaining summary statistics for the merged BAM of SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/samtools stats /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bam > /scratch/hp7d/lodestone-testing/SRR1582856/logs/7_bam_statistics/SRR1582856.samtools_stats_output.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/samtools flagstat /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bam > /scratch/hp7d/lodestone-testing/SRR1582856/logs/7_bam_statistics/SRR1582856.samtools_flagstat_output.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. calling variants from SRR1582856 using mpileup...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/core/anaconda2/4.4.0/bin/bcftools mpileup -Ou -f /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bam 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt | /apps/software/standard/core/anaconda2/4.4.0/bin/bcftools call --threads 20 --ploidy 1 -mv -Ov -o /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. regularising VCF for SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/project/som_infc_dis_mathers/tools/vcflib/bin/vcfallelicprimitives /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.regularised.vcf 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm --interactive=never /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.regularised.vcf /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. sorting, compressing and indexing VCF for SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
/apps/software/standard/compiler/gcc/7.1.0/vcftools/0.1.15/bin/vcf-sort /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.vcf 2>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mv /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.sorted.vcf /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf
/apps/software/standard/core/anaconda2/4.4.0/bin/bgzip /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf
/apps/software/standard/core/anaconda2/4.4.0/bin/tabix -p vcf /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. creating de novo assembly for SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
python /project/som_infc_dis_mathers/tools/SPAdes-3.11.1-Linux/bin/spades.py --cov-cutoff auto --careful -1 ./S1.fq.gz -2 ./S2.fq.gz -s /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz -t 20 -o /scratch/hp7d/lodestone-testing/SRR1582856/spades &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mv /scratch/hp7d/lodestone-testing/SRR1582856/spades/scaffolds.fasta /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa
rm -r /scratch/hp7d/lodestone-testing/SRR1582856/spades
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. creating summary statistics for the de novo assembly of SRR1582856...
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
mkdir -p /scratch/hp7d/lodestone-testing/SRR1582856/logs/8_assembly_statistics
python /project/som_infc_dis_mathers/tools/quast-5.0.1/bin/quast.py -r /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa -g /scratch/hp7d/lodestone-testing/SRR1582856/ref.gff -t 20 --circos --gene-finding --rna-finding --conserved-genes-finding --pe1 ./S1.fq.gz --pe2 ./S2.fq.gz --single /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz -o /scratch/hp7d/lodestone-testing/SRR1582856/logs/8_assembly_statistics /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa &>> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
rm ./S1.fq.gz ./S2.fq.gz /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.se.fq.gz
rm /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.amb /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.ann /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.bwt /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.fai /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.pac /scratch/hp7d/lodestone-testing/SRR1582856/ref.fa.sa /scratch/hp7d/lodestone-testing/SRR1582856/ref.gff
rm /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bam /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.bai
md5sum /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.fa.md5sum
md5sum /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf.gz > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf.gz.md5sum
md5sum /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf.gz.tbi > /scratch/hp7d/lodestone-testing/SRR1582856/SRR1582856.vcf.gz.tbi.md5sum
duration=$SECONDS
echo "--> total time taken to process SRR1582856: $(($duration / 60)) minutes and $(($duration % 60)) seconds
" >> /scratch/hp7d/lodestone-testing/log.SRR1582856.txt
exit 0
