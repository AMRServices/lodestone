#!/bin/bash

#########
### usage
#########
usage()
{
cat << EOF
usage: $0 [options]

This script generates lodestone submission commands
from demultiplexed output directory

Author: Hardik Parikh

OPTIONS:
	-h	Show this message
	-i	/FULL/PATH/TO/DEMUXDIR/
	-o	/FULL/PATH/TO/OUTPUTDIR/
EOF
}

if [ $# -eq 0 ]; then
        usage
        exit
fi


##################
### Accept options
##################
while getopts “hi:o:?” OPTION
do
        case $OPTION in
                h)      usage
                        exit 1;;
                i)      INDIR=$OPTARG;;
                o)      OUTDIR=$OPTARG;;
                ?)      usage
                        exit 1;;
        esac
done


########################
### List of commands ###
########################
RUNID=$(basename $INDIR)
touch ${RUNID}.lodestone_cmd.txt
touch ${RUNID}.lodestone_metrics_cmd.txt

for r1 in `find ${INDIR} -name "*_R1_*fastq.gz"`; 
do
	r2=${r1/_R1_/_R2_}
	r1file=$(basename $r1)
	cav=`echo $r1file | awk -F '_S' '{print $1}'`

	# echo the commands
	echo -e "lodestone.sh -c ${cav} -1 ${r1} -2 ${r2} -o ${OUTDIR}/${cav}" >> ${RUNID}.lodestone_cmd.txt 
	echo -e "collate_lodestone_metrics.sh -r ${RUNID} -c ${cav} -l ${OUTDIR}/${cav}" >> ${RUNID}.lodestone_metrics_cmd.txt 
done

### remove Undetermined from the list
grep -v "Undetermined" ${RUNID}.lodestone_cmd.txt > tmp
mv -f tmp ${RUNID}.lodestone_cmd.txt
grep -v "Undetermined" ${RUNID}.lodestone_metrics_cmd.txt > tmp
mv -f tmp ${RUNID}.lodestone_metrics_cmd.txt

echo "Done!"
exit 1
