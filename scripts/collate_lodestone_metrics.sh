#!/bin/bash


usage()
{
cat << EOF
usage: sh $0 options

This script collects various statistics from lodestone output folder 

Author: Hardik Parikh
Version: 1.0

OPTIONS:
	-h	Show this message
	-r	RUNID
	-c	CAVID
	-l	/FULL/PATH/TO/LODESTONE_OUTPUT/
EOF
}

if [ $# -eq 0 ]; then
        usage
        exit
fi


while getopts “hr:c:l:?” OPTION
do
        case $OPTION in
                h)      usage
                        exit 1;;
                r)      RUNID=$OPTARG;;
                c)      CAVID=$OPTARG;;
                l)      L_OUTDIR=$OPTARG;;
                ?)      usage
                        exit 1;;

        esac
done


### run fastqstats on raw r1/r2
### run fastqstats on hq r1/r2/se
DEMUX_DIR="/project/amr_services/demux"
#DEMUX_DIR="/scratch/hp7d/MathersLab/amr_lodestone_0120/demux"
mkdir -p ${L_OUTDIR}/logs/10_fastqstats
zcat ${DEMUX_DIR}/${RUNID}/${CAVID}*_R1_*.fastq.gz | /project/amr_services/scripts/fastqstats > ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R1.fastqstats.txt
zcat ${DEMUX_DIR}/${RUNID}/${CAVID}*_R2_*.fastq.gz | /project/amr_services/scripts/fastqstats > ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R2.fastqstats.txt
zcat ${L_OUTDIR}/${CAVID}_R1.fq.gz | /project/amr_services/scripts/fastqstats > ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R1.fastqstats.txt
zcat ${L_OUTDIR}/${CAVID}_R2.fq.gz | /project/amr_services/scripts/fastqstats > ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R2.fastqstats.txt
zcat ${L_OUTDIR}/${CAVID}.se.fq.gz | /project/amr_services/scripts/fastqstats > ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_se.fastqstats.txt

# collect read counts
rawreads_for=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R1.fastqstats.txt | grep "Number of reads" | awk '{print $NF}'`
rawreads_rev=${rawreads_for}
hqreads_for=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R1.fastqstats.txt | grep "Number of reads" | awk '{print $NF}'`
hqreads_rev=${hqreads_for} 
hqreads_se=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_se.fastqstats.txt | grep "Number of reads" | awk '{print $NF}'`
# collect read length
rawreads_len_for=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R1.fastqstats.txt | grep "Avg read length" | awk '{printf("%.2f\n", $NF)}'`
rawreads_len_rev=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R2.fastqstats.txt | grep "Avg read length" | awk '{printf("%.2f\n", $NF)}'`
hqreads_len_for=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R1.fastqstats.txt | grep "Avg read length" | awk '{printf("%.2f\n", $NF)}'`
hqreads_len_rev=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R2.fastqstats.txt | grep "Avg read length" | awk '{printf("%.2f\n", $NF)}'`
hqreads_len_se=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_se.fastqstats.txt | grep "Avg read length" | awk '{printf("%.2f\n", $NF)}'`
# collect avgQ
rawreads_qual_for=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R1.fastqstats.txt | grep "Avg read quality" | awk '{printf("%.2f\n", $NF)}'`
rawreads_qual_rev=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_raw_R2.fastqstats.txt | grep "Avg read quality" | awk '{printf("%.2f\n", $NF)}'`
hqreads_qual_for=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R1.fastqstats.txt | grep "Avg read quality" | awk '{printf("%.2f\n", $NF)}'`
hqreads_qual_rev=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_R2.fastqstats.txt | grep "Avg read quality" | awk '{printf("%.2f\n", $NF)}'`
hqreads_qual_se=`cat ${L_OUTDIR}/logs/10_fastqstats/${CAVID}_hq_se.fastqstats.txt | grep "Avg read quality" | awk '{printf("%.2f\n", $NF)}'`

# kaiju hits
top_sp_pe=$(awk -F '\t' 'NR==3 {print $3}' ${L_OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
top_sp_perc_pe=$(awk -F '\t' 'NR==3 {print $1}' ${L_OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
top_sp_perc_pe=`echo $top_sp_perc_pe | awk '{printf("%.2f", $1)}'`
sec_sp_pe=$(awk -F '\t' 'NR==4 {print $3}' ${L_OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
sec_sp_perc_pe=$(awk -F '\t' 'NR==4 {print $1}' ${L_OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_species.pe.txt)
sec_sp_perc_pe=`echo $sec_sp_perc_pe | awk '{printf("%.2f", $1)}'`

# candidate taxon
taxlevel=`grep "taxon ID" ${L_OUTDIR}/log.${CAVID}.txt | awk '{print $10}'`
if [[ $taxlevel == "species" ]]; then
    taxaname=${top_sp_pe}
else
    taxaname=$(awk -F '\t' 'NR==3 {print $3}' ${L_OUTDIR}/logs/5_kaiju_report/${CAVID}.kaiju_genus.pe.txt)
fi
taxid=`grep "taxon ID" ${L_OUTDIR}/log.${CAVID}.txt | awk '{print $13}'`
num_genomes=`grep "taxon ID" ${L_OUTDIR}/log.${CAVID}.txt | awk '{print $6}'`

# mash results
mash_closest_ref_acc=`grep "NCBI_Assembly" ${L_OUTDIR}/ref.gff | awk -F ':' '{print $2}'`
mash_closest_ref_str=`grep ${mash_closest_ref_acc} ${L_OUTDIR}/logs/6_ref_genomes_downloaded/assembly_summary.txt | awk -F '\t' '{print $8,$9}'`

# chrom seq id
chrom=`head -1 ${L_OUTDIR}/ref.fa.fai | awk -F '\t' '{print $1}'`
chromlen=`head -1 ${L_OUTDIR}/ref.fa.fai | awk -F '\t' '{print $2}'`

# bwa results
perc_mapped=`grep "mapped (" ${L_OUTDIR}/logs/8_bam_statistics/*flagstat_output.txt | awk -F '(' '{print $2}' | awk '{print $1}'`

# genome coverage
genomeCoverageBed -d -ibam ${L_OUTDIR}/${CAVID}.bam -g ${L_OUTDIR}/ref.fa.fai > ${L_OUTDIR}/logs/8_bam_statistics/${CAVID}.bedtools.genomecov.txt
nonzero_pos=`cat ${L_OUTDIR}/logs/8_bam_statistics/${CAVID}.bedtools.genomecov.txt | grep ${chrom} | awk -F '\t' '$3 != 0' | wc -l | awk '{print $1}'`
perc_genome_covered=`awk -v nz=${nonzero_pos} -v len=${chromlen} 'BEGIN {printf("%.2f", nz*100.0/len)}'`
coverage=`cat ${L_OUTDIR}/logs/8_bam_statistics/${CAVID}.bedtools.genomecov.txt | grep ${chrom} | awk -F '\t' '$3 != 0' | awk -F '\t' '{s+=$3}END{printf("%.2f", s/NR)}'`

# VCF stats
snps=`zcat ${L_OUTDIR}/${CAVID}.vcf.gz | grep -v "^#" | grep -v "TYPE" | grep -v "INDEL"| grep ${chrom} | wc -l`
hqsnps=`zcat ${L_OUTDIR}/${CAVID}.vcf.gz | /project/som_infc_dis_mathers/tools/bcftools/misc/vcfutils.pl varFilter -Q 30 -d 20 | grep -v "^#" | grep -v "TYPE" | grep -v "INDEL"| grep ${chrom} | wc -l`

# Quast stats
num_contigs=`grep "^# contigs" ${L_OUTDIR}/logs/9_assembly_statistics/quast_out/report.txt| grep -v "(" | awk '{print $NF}'`
num_contigs_1kb=`grep "contigs (>= 1000 bp)" ${L_OUTDIR}/logs/9_assembly_statistics/quast_out/report.txt | awk '{print $NF}'`
tot_len=`grep "Total length" ${L_OUTDIR}/logs/9_assembly_statistics/quast_out/report.txt | grep -v "(" | awk '{print $NF}'`
gc=`grep "GC (%)" ${L_OUTDIR}/logs/9_assembly_statistics/quast_out/report.txt | awk '{print $NF}'`
n50=`grep "N50" ${L_OUTDIR}/logs/9_assembly_statistics/quast_out/report.txt | awk '{print $NF}'`

# Lodestone Erros/Warnings
if [[ -f ${L_OUTDIR}/${CAVID}.vcf.gz.md5sum || -f ${L_OUTDIR}/logs/9_assembly_statistics/quast_out/report.txt ]]; then
	err="Pipeline finished."
else
	err="ERRORS. Please check ${L_OUTDIR}/log.${CAVID}.txt!"
fi

stime=`grep "^--" ${L_OUTDIR}/log.${CAVID}.txt | head -1 | cut -d "." -f 1 | awk -F ">" '{print $2}'`
etime=`grep "^--" ${L_OUTDIR}/log.${CAVID}.txt | tail -1 | cut -d "." -f 1 | awk -F ">" '{print $2}'`
proctime=`echo $(( $(date '+%s' --date="${etime}") - $(date '+%s' --date="${stime}") )) | awk '{print int($1/60)":"int($1%60)}'`

# print output
OUTFILE=${L_OUTDIR}/${CAVID}_lodestone_summary.tsv
touch ${OUTFILE}
#RUNID=`dirname ${L_OUTDIR} | awk -F '/' '{print $NF}'`

echo -e "RunID\tSampleID\tRawReads\tHQReads\tRawReadLen(Avg)\tHQReadLen(Avg)\tRawReadQual(Avg)\tHQReadQual(Avg)\tKaiju_PE_TopHit\tKaiju_PE_SecHit\tCandidateTaxon\tCandidateTaxon_NCBITaxID\tNumGenomeDwl\tMASH_ClosestRef\tMASH_ClosestRef_AccID\tBWA_ReadsMapped(%)\tGenomeBasesCovered\tGenomeCoverage(%)\tDepthofCoverage(X)\tSNPs\tHQ_SNPs(QUAL>=30, DP>=20)\tNumContigs\tNumContigs(>=1kb)\tTotalLength\tGC(%)\tN50\tLodestoneError\tLodestone_TotalProcTime(min:sec)" > ${OUTFILE}

echo -e "${RUNID}\t${CAVID}\t${rawreads_for}|${rawreads_rev}\t${hqreads_for}|${hqreads_rev}|${hqreads_se}\t${rawreads_len_for}|${rawreads_len_rev}\t${hqreads_len_for}|${hqreads_len_rev}|${hqreads_len_se}\t${rawreads_qual_for}|${rawreads_qual_rev}\t${hqreads_qual_for}|${hqreads_qual_rev}|${hqreads_qual_se}\t${top_sp_perc_pe}% ${top_sp_pe}\t${sec_sp_perc_pe}% ${sec_sp_pe}\t${taxname}\t${taxid}\t${num_genomes}\t${mash_closest_ref_str}\t${mash_closest_ref_acc}\t${perc_mapped}\t${nonzero_pos}\t${perc_genome_covered}%\t${coverage}\t${snps}\t${hqsnps}\t${num_contigs}\t${num_contigs_1kb}\t${tot_len}\t${gc}\t${n50}\t${err}\t${proctime}" >> ${OUTFILE}
