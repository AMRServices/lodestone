#!/bin/bash
PATH=$PATH:/mnt/microbio/HOMES/steveb/programs/circos-0.69-6/bin
PATH=$PATH:/mnt/microbio/HOMES/steveb/programs/barrnap/bin
PATH=$PATH:/mnt/microbio/HOMES/steveb/programs/samtools-1.7
PATH=$PATH:/mnt/microbio/HOMES/steveb/programs/bcftools
PATH=$PATH:/mnt/microbio/HOMES/steveb/programs/ncbi-blast-2.2.25+/bin
PATH=$PATH:/mnt/microbio/HOMES/steveb/programs/bwa-0.7.17
SECONDS=0
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "-->$current_date_time. beginning analysis of SRR974838
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "
-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. testing software availability before beginning..." >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Aspera ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/home/ndm.local/steveb/.aspera/connect/bin/ascp --version >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### FastQC ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/FastQC/fastqc --version >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### cutadapt ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/home/ndm.local/steveb/.local/bin/cutadapt --version >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### TrimGalore ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/TrimGalore-0.5.0/trim_galore --version >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Trimmomatic ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/Trimmomatic-0.38/trimmomatic-0.38.jar -version >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Centrifuge ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/centrifuge/centrifuge -version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Kaiju ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/kaiju/bin/kaiju -h &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### McCortex ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/mccortex/bin/mccortex31 -h &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Atlas ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/home/ndm.local/steveb/.local/bin/atlas --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Mykrobe ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/home/ndm.local/steveb/.local/bin/mykrobe --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Mash ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/mash-Linux64-v2.1/mash --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### seqtk ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/seqtk/seqtk &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Entrez Direct ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/edirect/efetch --help &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### BWA ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bwa-0.7.17/bwa &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### sambamba ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/sambamba --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### SAMtools ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/samtools-1.7/samtools &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### BCFtools ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bcftools/bcftools --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Picard Tools ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar MarkDuplicates --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### VCFlib ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/vcflib/bin/vcfoverlay --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### VCFtools ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bin/vcf-sort --help &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### bgzip ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bin/bgzip --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### tabix ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bin/tabix --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### SPAdes ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
python /mnt/microbio/HOMES/steveb/programs/SPAdes-3.13.0-Linux/bin/spades.py --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Quast ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
python /mnt/microbio/HOMES/steveb/programs/quast-5.0.1/quast.py --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Circos ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/circos-0.69-6/bin/circos --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Prokka ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/prokka/bin/prokka --version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### TETyper ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/home/ndm.local/steveb/.local/bin/TETyper.py -h &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### BLAST+ ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/ncbi-blast-2.2.25+/bin/blastn -version &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "
-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. confirming database availability..." >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Database employed by Centrifuge ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
date_last_modified="`date -r /mnt/microbio/HOMES/steveb/programs/centrifuge/db.human_and_hiv`"
echo "/mnt/microbio/HOMES/steveb/programs/centrifuge/db.human_and_hiv. last modified $date_last_modified" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
echo "
### Databases employed by Kaiju ###
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
date_last_modified="`date -r /mnt/microbio/HOMES/steveb/programs/kaiju/db.p/kaiju_db.fmi`"
echo "/mnt/microbio/HOMES/steveb/programs/kaiju/db.p/kaiju_db.fmi. last modified $date_last_modified" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
date_last_modified="`date -r /mnt/microbio/HOMES/steveb/programs/kaiju/db.e/kaiju_db_nr_euk.fmi`"
echo "/mnt/microbio/HOMES/steveb/programs/kaiju/db.e/kaiju_db_nr_euk.fmi. last modified $date_last_modified" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "
-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. beginning pipeline...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "the last command run - ${last_command} - finished with exit code $?."' EXIT
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838
cd /mnt/microbio/HOMES/steveb/sp3_output/SRR974838
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. downloading fastqs for SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
if curl --head --fail --silent "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR974/SRR974838/SRR974838_1.fastq.gz" >/dev/null; then wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR974/SRR974838/SRR974838_1.fastq.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log; else echo "ERROR: unable to find file at URL ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR974/SRR974838/SRR974838_1.fastq.gz; is this because SRR974838 does not contain PE fqs?" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
if curl --head --fail --silent "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR974/SRR974838/SRR974838_2.fastq.gz" >/dev/null; then wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR974/SRR974838/SRR974838_2.fastq.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log; else echo "ERROR: unable to find file at URL ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR974/SRR974838/SRR974838_2.fastq.gz; is this because SRR974838 does not contain PE fqs?" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
if [ ! -e $fq_1 ]; then echo "ERROR: expected PE reads but unable to find /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
if [ ! -e $fq_2 ]; then echo "ERROR: expected PE reads but unable to find /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
if curl --head --fail --silent "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz" >/dev/null; then wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz -O /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log; else echo "ERROR: a reference genome was specified but we were unable to find the file at URL ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
if curl --head --fail --silent "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz" >/dev/null; then wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz -O /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.gff.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log; else echo "ERROR: a reference GFF was specified but we were unable to find the file at URL ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
if curl --head --fail --silent "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_protein.faa.gz" >/dev/null; then wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_protein.faa.gz -O /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.faa.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log; else echo "ERROR: a reference FAA was specified but we were unable to find the file at URL ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_protein.faa.gz" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log && exit 1; fi
gunzip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.gz
gunzip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.gff.gz
gunzip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.faa.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running initial FastQC on SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/FastQC/fastqc /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/1_initial_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1_fastqc.html /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/1_initial_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2_fastqc.html /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/1_initial_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1_fastqc.zip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/1_initial_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2_fastqc.zip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/1_initial_fastqc_report
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running cutadapt/TrimGalore on SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/TrimGalore-0.5.0/trim_galore --path_to_cutadapt /home/ndm.local/steveb/.local/bin/cutadapt --stringency 5 --trim-n --max_n 5 --paired --length 0 -o /mnt/microbio/HOMES/steveb/sp3_output/SRR974838 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/2_trimgalore_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz_trimming_report.txt /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/2_trimgalore_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz_trimming_report.txt /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/2_trimgalore_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1_val_1.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2_val_2.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/3_trimmomatic_report
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running Trimmomatic on the adapter-trimmed reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -threads 6 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.forward_paired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.forward_unpaired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.reverse_paired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.reverse_unpaired.fq.gz TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100 >/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/3_trimmomatic_report/SRR974838.trimmomatic_log.txt 2>&1
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.forward_paired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.reverse_paired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz
cat /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.forward_unpaired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.reverse_unpaired.fq.gz > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.forward_unpaired.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.reverse_unpaired.fq.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running FastQC on the cleaned and adapter-trimmed reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/FastQC/fastqc /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1_fastqc.html /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2_fastqc.html /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1_fastqc.zip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2_fastqc.zip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se_fastqc.html /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se_fastqc.zip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/4_final_fastqc_report
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. running Centrifuge on the cleaned and adapter-trimmed reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/centrifuge/centrifuge -p 6 -x /mnt/microbio/HOMES/steveb/programs/centrifuge/db.human_and_hiv/db_human_and_hiv -1 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz -2 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz -U /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz --report-file /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads.txt -S /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.centrifuge_read_classification &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
awk '$2=="unclassified" { print $1 }' /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.centrifuge_read_classification > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.non_human_read_list
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.centrifuge_read_classification
awk 'NR==1 { print $0 }' /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads.txt > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads_excl_hiv.txt
grep "Homo sapiens" /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads.txt >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads_excl_hiv.txt
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads_excl_hiv.txt /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.number_of_human_reads.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing reads classified as human from the cleaned and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/seqtk/seqtk subseq /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.non_human_read_list | gzip > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.excl_human.fastq.gz
/mnt/microbio/HOMES/steveb/programs/seqtk/seqtk subseq /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.non_human_read_list | gzip > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.excl_human.fastq.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing reads classified as human from the cleaned and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/seqtk/seqtk subseq /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.non_human_read_list > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.excl_human.fastq.gz
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/5_centrifuge_report_on_number_of_human_reads/SRR974838.non_human_read_list
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.excl_human.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.excl_human.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.excl_human.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. indexing reference genome...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bwa-0.7.17/bwa index /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. aligning the quality- and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bwa-0.7.17/bwa mem -R '@RG\tID:group_SRR974838\tSM:sample_SRR974838\tPL:Illumina\tLIB:lib_SRR974838\tPU:unit_SRR974838' -t 6 -M /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log | /mnt/microbio/HOMES/steveb/programs/sambamba view -S -f bam -t 6 -o /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam /dev/stdin
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. cleaning BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar CleanSam INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.cleaned.paired-end.bam TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.cleaned.paired-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. fixing mate information in the BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar FixMateInformation INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.fixmated.paired-end.bam TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.fixmated.paired-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. sorting BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar SortSam INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.paired-end.bam SORT_ORDER=coordinate TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. duplicate-marking BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar MarkDuplicates INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.paired-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bam METRICS_FILE=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.metrics ASSUME_SORTED=true TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing supplementary and non-primary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed paired-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/samtools-1.7/samtools view -b -F 2048 -F 256 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bam > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.rm_supplementary.bam 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.rm_supplementary.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bam
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar BuildBamIndex INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bam &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.paired-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.paired-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.metrics
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. aligning the quality- and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bwa-0.7.17/bwa mem -R '@RG\tID:group_SRR974838\tSM:sample_SRR974838\tPL:Illumina\tLIB:lib_SRR974838\tPU:unit_SRR974838' -t 6 -M /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log | /mnt/microbio/HOMES/steveb/programs/sambamba view -S -f bam -t 6 -o /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam /dev/stdin
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. cleaning BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar CleanSam INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.cleaned.single-end.bam TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.cleaned.single-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. fixing mate information in the BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar FixMateInformation INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.fixmated.single-end.bam TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.fixmated.single-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. sorting BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar SortSam INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.single-end.bam SORT_ORDER=coordinate TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. duplicate-marking BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar MarkDuplicates INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.single-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bam METRICS_FILE=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.metrics ASSUME_SORTED=true TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. removing supplementary and non-primary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed single-end reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/samtools-1.7/samtools view -b -F 2048 -F 256 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bam > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.rm_supplementary.bam 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.rm_supplementary.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bam
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar BuildBamIndex INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bam &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.unsorted.single-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.single-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.metrics
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. merging PE and SE BAMs for the aligned, classified, quality- and adapter-trimmed reads from SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar MergeSamFiles SORT_ORDER='coordinate' INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bam INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bam OUTPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bam TMP_DIR=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838 &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
java -jar /mnt/microbio/HOMES/steveb/programs/picard.jar BuildBamIndex INPUT=/mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bam &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.paired-end.bai /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.single-end.bai
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/6_bam_statistics
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. obtaining summary statistics for the merged BAM of SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/samtools-1.7/samtools stats /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bam > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/6_bam_statistics/SRR974838.samtools_stats_output.txt
/mnt/microbio/HOMES/steveb/programs/samtools-1.7/samtools flagstat /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bam > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/6_bam_statistics/SRR974838.samtools_flagstat_output.txt
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. calling variants from SRR974838 using mpileup...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bcftools/bcftools mpileup -Ou -f /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bam 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log | /mnt/microbio/HOMES/steveb/programs/bcftools/bcftools call --threads 6 --ploidy 1 -mv -Ov -o /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. regularising VCF for SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/vcflib/bin/vcfallelicprimitives /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.regularised.vcf 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm --interactive=never /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.regularised.vcf /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. sorting, compressing and indexing VCF for SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
/mnt/microbio/HOMES/steveb/programs/bin/vcf-sort /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.vcf 2>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.sorted.vcf /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf
/mnt/microbio/HOMES/steveb/programs/bin/bgzip /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf
/mnt/microbio/HOMES/steveb/programs/bin/tabix -p vcf /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf.gz
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. creating de novo assembly for SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
python /mnt/microbio/HOMES/steveb/programs/SPAdes-3.13.0-Linux/bin/spades.py --cov-cutoff auto --careful -1 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz -2 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz -s /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz -t 6 -o /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/spades &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mv /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/spades/scaffolds.fasta /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.fa
rm -r /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/spades
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. creating summary statistics for the de novo assembly of SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
mkdir -p /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/7_quast_report_on_de_novo_assembly
python /mnt/microbio/HOMES/steveb/programs/quast-5.0.1/quast.py -r /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa -g /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.gff -t 6 --circos --gene-finding --rna-finding --conserved-genes-finding --pe1 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz --pe2 /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz --single /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz -o /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/7_quast_report_on_de_novo_assembly /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.fa &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
elapsed_time=$SECONDS
echo "-->$current_date_time. elapsed time: $(($elapsed_time / 60)) minutes and $(($elapsed_time % 60)) seconds. annotating the de novo assembly of SRR974838...
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
species_name_for_prokka=$( echo $top_species_hit_pe | perl -pe 's/(.*?) (.*?)/$2/' )
strain_name_for_prokka=$( echo $infraspecific_name | perl -pe 's/strain\=//' )
/mnt/microbio/HOMES/steveb/programs/prokka/bin/prokka --cpus 6 --kingdom Bacteria --outdir /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/logs/8_prokka_annotation_of_de_novo_assembly --prefix SRR974838 --locustag SRR974838 --compliant --genus "$top_genus_hit_pe" --species "$species_name_for_prokka" --strain "$strain_name_for_prokka" --proteins /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.faa --usegenus --rfam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.fa &>> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_1.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838_2.fastq.gz /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.se.fq.gz
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.amb /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.ann /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.bwt /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.fai /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.pac /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.fa.sa /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/ref.gff
rm /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bam /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.bai
md5sum /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.fa > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.fa.md5sum
md5sum /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf.gz > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf.gz.md5sum
md5sum /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf.gz.tbi > /mnt/microbio/HOMES/steveb/sp3_output/SRR974838/SRR974838.vcf.gz.tbi.md5sum
duration=$SECONDS
echo "--> total time taken to process SRR974838: $(($duration / 60)) minutes and $(($duration % 60)) seconds
" >> /mnt/microbio/HOMES/steveb/sp3_output/SRR974838.log
exit 0
