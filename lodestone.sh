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
Version: 0.2

OPTIONS:
        -h      Show this message
        -c      CAVID
        -1      FULL/PATH/TO/ForwardRawReads (*_R1.fq.gz)
        -2      FULL/PATH/TO/ReverseRawReads (*_R2.fq.gz)
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
while getopts “hc:1:2:o:?” OPTION
do
        case $OPTION in
                h)      usage
                        exit 1;;
                c)      CAVID=$OPTARG;;
                1)      IN_R1FILE=$OPTARG;;
                2)      IN_R2FILE=$OPTARG;;
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

SECONDS=0
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. beginning analysis of ${CAVID}\n" >> ${LOGFILE}
cp ${IN_R1FILE} ${R1FILE}
cp ${IN_R2FILE} ${R2FILE}


################
### Dependencies
################
# fastqc, trim_galore, trimmomatic
# kaiju, mash
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


#####################
### Speciation: Kaiju 
#####################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running Kaiju on high-quality paired-end reads from ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/5_kaiju_report
kaiju -z ${THREADS} -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -f /project/som_infc_dis_mathers/tools/kaiju/kaijudb/kaiju_db.fmi -i ${R1FILE} -j ${R2FILE} -a greedy -e 5 -E 0.05 -o ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_report.pe.txt
kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_report.pe.txt -r genus -o ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_genus.pe.txt
kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_report.pe.txt -r species -o ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running Kaiju on high-quality single-end reads from ${CAVID}...\n" >> ${LOGFILE}
kaiju -z ${THREADS} -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -f /project/som_infc_dis_mathers/tools/kaiju/kaijudb/kaiju_db.fmi -i ${OUTDIR}/${CAVID}.se.fq.gz -a greedy -e 5 -E 0.05 -o ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_report.se.txt
kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_report.se.txt -r genus -o ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_genus.se.txt
kaijuReport -t /project/som_infc_dis_mathers/tools/kaiju/kaijudb/nodes.dmp -n /project/som_infc_dis_mathers/tools/kaiju/kaijudb/names.dmp -i ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_report.se.txt -r species -o ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.se.txt


########################
### Download Ref Genomes
########################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. downloading reference genomes for top species/genus identified by Kaiju for ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/6_ref_genomes_downloaded

## Check if top-hit for PE and SE match, if not - print warning, proceed with top-hit from PE
top_genus_pe_hit=$(awk -F '\t' 'NR==3 {print $3}' ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_genus.pe.txt)
top_sp_pe_hit=$(awk -F '\t' 'NR==3 {print $3}' ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
top_sp_se_hit=$(awk -F '\t' 'NR==3 {print $3}' ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.se.txt)
if [[ "$top_sp_pe_hit" != "$top_sp_se_hit" ]]; then 
	echo -e "\n-->$current_date_time. WARNING: candidate species for ${CAVID} as predicted from PE and SE reads DO NOT match; proceeding with paired-end hit ...\n" >> ${LOGFILE}
else
	echo -e "\n-->$current_date_time. candidate species for ${CAVID} as predicted from PE and SE reads MATCH!\n" >> ${LOGFILE}
fi

## Pick the taxonomic level to download references for 
## If top species has <50% reads classified, AND second hit is from same genus - 
## then move to genus-level 
tax_level_to_dwl="species"
top_sp_pe_hit=$(awk -F '\t' 'NR==3 {print $3}' ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
top_sp_reads=$(awk -F '\t' 'NR==3 {print $1}' ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
top_hit_genus=`echo $top_sp_pe_hit | awk '{print $1}'`
sec_sp_pe_hit=$(awk -F '\t' 'NR==4 {print $3}' ${OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
sec_hit_genus=`echo $sec_sp_pe_hit | awk '{print $1}'`
if ([[ "$top_hit_genus" == "$sec_hit_genus" ]] && [[ "$top_sp_reads" < 50.0 ]]); then
	tax_level_to_dwl="genus"
fi

## Get TaxID for top hit
SPECIES_TAXID=$(esearch -db taxonomy -query "$top_sp_pe_hit" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId)
GENUS_TAXID=$(esearch -db taxonomy -query "$top_genus_pe_hit" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId)

## Download genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt &>> ${LOGFILE}
date_last_modified="`date -r assembly_summary.txt`"
echo -e "\n\n--> downloaded ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt. last modified $date_last_modified" >> ${LOGFILE}

if [[ "$tax_level_to_dwl" == "species" ]]; then
	awk -v sptaxid=${SPECIES_TAXID} -F "\t" '$7==sptaxid && $12=="Complete Genome" && $11=="latest"{print $20}' ${OUTDIR}/assembly_summary.txt > ${OUTDIR}/ftpdirpaths.txt
	awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ${OUTDIR}/ftpdirpaths.txt > ${OUTDIR}/ftpfilepaths.txt
	number_of_complete_genomes="`wc -l < ${OUTDIR}/ftpfilepaths.txt`"
	echo -e "\n\n--> assembly_summary.txt parsed to identify $number_of_complete_genomes complete genomes for species taxon ID $SPECIES_TAXID" >> ${LOGFILE}
	if [[ "$number_of_complete_genomes" == 0 ]]; then
		awk -F '\t' -v genus=${top_genus_pe_hit} '$8 ~ genus' ${OUTDIR}/assembly_summary.txt | awk -F '\t' '$12=="Complete Genome" && $11=="latest" {print $20}' > ${OUTDIR}/ftpdirpaths.txt
		awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ${OUTDIR}/ftpdirpaths.txt > ${OUTDIR}/ftpfilepaths.txt
		number_of_complete_genomes="`wc -l < ${OUTDIR}/ftpfilepaths.txt`"
		echo -e "\n\n--> No reference genome for species taxon ID ${SPECIES_TAXID}, expanding the search to genus-level..." >> ${LOGFILE}
		echo -e "--> assembly_summary.txt parsed to identify $number_of_complete_genomes complete genomes for genus taxon ID $GENUS_TAXID" >> ${LOGFILE}
	fi
else
	awk -F '\t' -v genus=${top_genus_pe_hit} '$8 ~ genus' ${OUTDIR}/assembly_summary.txt | awk -F '\t' '$12=="Complete Genome" && $11=="latest" {print $20}' > ${OUTDIR}/ftpdirpaths.txt
	awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ${OUTDIR}/ftpdirpaths.txt > ${OUTDIR}/ftpfilepaths.txt
	number_of_complete_genomes="`wc -l < ${OUTDIR}/ftpfilepaths.txt`"
	echo -e "\n\n--> Kaiju PE top-hit < 50% reads AND second-hit from different genus: expanding the search to genus-level..." >> ${LOGFILE}
	echo -e "--> assembly_summary.txt parsed to identify $number_of_complete_genomes complete genomes for genus taxon ID $GENUS_TAXID" >> ${LOGFILE}
fi

## Exit pipeline if no reference genomes
if [[ "$number_of_complete_genomes" == 0 ]]; then
	current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
	echo -e "\n\n-->$current_date_time. No Reference Genomes Downloaded for ${CAVID}. Unable to proceed!\n"
	exit 1
fi

echo "NUMBER OF COMPLETE GENOMES OBTAINED AFTER PARSING ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt (DOWNLOADED $date_last_modified): $number_of_complete_genomes" >> ${OUTDIR}/logs/6_ref_genomes_downloaded/${CAVID}.reference_genome_download_log.txt
echo "GENOME URLS:" >> ${OUTDIR}/logs/6_ref_genomes_downloaded/${CAVID}.reference_genome_download_log.txt
cat ${OUTDIR}/ftpfilepaths.txt >> ${OUTDIR}/logs/6_ref_genomes_downloaded/${CAVID}.reference_genome_download_log.txt
mkdir -p ${OUTDIR}/temp_dir_genome_storage
wget -i ${OUTDIR}/ftpfilepaths.txt -P ${OUTDIR}/temp_dir_genome_storage &>> ${LOGFILE}


########
### MASH
########
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. running MASH to identify closest reference genome for ${CAVID}...\n" >> ${LOGFILE}
mkdir -p ${OUTDIR}/logs/7_mash_report
find ${OUTDIR}/temp_dir_genome_storage -maxdepth 1 -type f > ${OUTDIR}/fa_filepaths_for_mash_sketch.txt
cat ${R1FILE} ${R2FILE} ${OUTDIR}/${CAVID}.se.fq.gz > ${OUTDIR}/${CAVID}.concatenated_reads.fq.gz
mash sketch -p ${THREADS} -l ${OUTDIR}/fa_filepaths_for_mash_sketch.txt -o ${OUTDIR}/${CAVID}.fa &>> ${LOGFILE}
mash sketch -p ${THREADS} -o ${OUTDIR}/${CAVID}.fq ${OUTDIR}/${CAVID}.concatenated_reads.fq.gz &>> ${LOGFILE}
rm ${OUTDIR}/${CAVID}.concatenated_reads.fq.gz
mash dist -p 20 ${OUTDIR}/${CAVID}.fq.msh ${OUTDIR}/${CAVID}.fa.msh > ${OUTDIR}/logs/7_mash_report/${CAVID}.mash-dist_output.tsv
closest_mash_hit=$(basename $(sort -g -k3 ${OUTDIR}/logs/7_mash_report/${CAVID}.mash-dist_output.tsv | head -1 | cut -f2))
closest_mash_hit_root=${closest_mash_hit%_genomic.fna.gz}
path_to_closest_mash_hit=$(grep $closest_mash_hit_root ${OUTDIR}/assembly_summary.txt | head -1 | awk -F '\t' '{print $20}')
FTPPATHG=$path_to_closest_mash_hit/$closest_mash_hit_root'_genomic.fna.gz'
FTPPATHGFF=$path_to_closest_mash_hit/$closest_mash_hit_root'_genomic.gff.gz'
wget $FTPPATHG -O ${OUTDIR}/ref.fa.gz &>> ${LOGFILE}
gunzip ${OUTDIR}/ref.fa.gz
wget $FTPPATHGFF -O ${OUTDIR}/ref.gff.gz &>> ${LOGFILE}
gunzip ${OUTDIR}/ref.gff.gz
echo "URL FOR CLOSEST MASH HIT: $FTPPATHG" >> ${OUTDIR}/logs/7_mash_report/${CAVID}.reference_genome_download_log.txt
mv ${OUTDIR}/assembly_summary.txt ${OUTDIR}/logs/6_ref_genomes_downloaded/
rm ${OUTDIR}/fa_filepaths_for_mash_sketch.txt
rm ${OUTDIR}/${CAVID}.fq.msh ${OUTDIR}/${CAVID}.fa.msh
rm -rf ${OUTDIR}/temp_dir_genome_storage


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
cp ${OUTDIR}/spades/contigs.fasta ${OUTDIR}/${CAVID}.fa
### check this ###
quast.py -t ${THREADS} --pe1 ${R1FILE} --pe2 ${R2FILE} --single ${OUTDIR}/${CAVID}.se.fq.gz -o ${OUTDIR}/logs/9_assembly_statistics/quast_out ${OUTDIR}/${CAVID}.fa &>> ${LOGFILE}


######################
### MultiQC Report ###
######################
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo -e "\n\n-->$current_date_time. Generating MultiQC report for ${CAVID}...\n" >> ${LOGFILE}
cd ${OUTDIR}
multiqc ${OUTDIR}


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
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
exit 0
