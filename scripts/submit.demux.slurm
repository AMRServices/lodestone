#!/bin/bash
#SBATCH --job-name=demux		# Job name
#SBATCH --cpus-per-task=20					# Number of CPU cores per task
#SBATCH --mem=60gb							# Job Memory
#SBATCH --time=12:00:00						# Time limit hrs:min:sec
#SBATCH --output=demux_%A.out				# Standard output log
#SBATCH --error=demux_%A.err				# Standard error log
#SBATCH -A amr_services						# allocation groups
#SBATCH -p standard							# slurm queue
pwd; hostname; date

### Load modules 
module purge
module use /project/som_infc_dis_mathers/modulefiles
module load bcl2fastq2
module load fastqc
module load multiqc

### command
 
sh  demuxMiSeq.sh 
	-i /project/amr_services/miseq/181105_M70720_0007_000000000-C574B \
	-o /scratch/mst3k/amr_services/demux/ \
	-t 20

date
