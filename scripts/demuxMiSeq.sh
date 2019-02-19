#!/bin/bash

#########
### usage
#########
usage()
{
cat << EOF
usage: sh $0 options

This script uses the bcl2fastq2 conversion software 
and run the BCL conversion and demultiplexing process 
for MiSeq raw data

The script also runs FastQC on each pair of reads
and summarizes the results using MultiQC

Author: Hardik Parikh
Version: bcl2fastq v2.20.0.422

OPTIONS:
	-h	Show this message
	-i	/FULL/PATH/MiSeq/RUNDATE/
	-o	/FULL/PATH/TO/OUTBASEDIR/
	-t	number of threads 
EOF
}

if [ $# -eq 0 ]; then
        usage
        exit
fi


##################
### Accept options
##################
while getopts “hi:o:t:?” OPTION
do
        case $OPTION in
                h)      usage
                        exit 1;;
                i)      INDIR=$OPTARG;;
                o)      OUTBASEDIR=$OPTARG;;
                t)      THREADS=$OPTARG;;
                ?)      usage
                        exit 1;;
        esac
done


###################
### Directory Names
###################
RUN_DIR_NAME=`echo $INDIR | awk -F '/' '{print $NF}'`
OUTDIR=${OUTBASEDIR}/${RUN_DIR_NAME}

##############
### bcl2fastq2
##############

bcl2fastq \
	--input-dir	${INDIR}/Data/Intensities/BaseCalls/ \
	--runfolder-dir ${INDIR} \
	--output-dir ${OUTDIR} \
	--sample-sheet ${INDIR}/SampleSheet.csv \
	--loading-threads ${THREADS} \
	--processing-threads ${THREADS}


##########
### FastQC
##########

mkdir ${OUTDIR}/FastQC
fastqc -t 20 -o ${OUTDIR}/FastQC ${OUTDIR}/*.fastq.gz

###########
### MultiQC
###########

cd ${OUTDIR}
multiqc .

###########
echo "Done!"
exit 1
