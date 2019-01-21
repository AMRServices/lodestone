#!/bin/bash
#SBATCH --job-name=lodestone_CAV001		# Job name
#SBATCH --cpus-per-task=20					# Number of CPU cores per task
#SBATCH --mem=80gb							# Job Memory
#SBATCH --time=12:00:00						# Time limit hrs:min:sec
#SBATCH --output=lodestone_CAV001_%A.out				# Standard output log
#SBATCH --error=lodestone_CAV001_%A.err				# Standard error log
#SBATCH -A amr_services						# allocation groups
#SBATCH -p standard							# slurm queue
pwd; hostname; date

### Load modules 
module purge
module use /project/som_infc_dis_mathers/modulefiles
module load lodestone


### command 
sh lodestone.sh \
	-c CAV001 \
	-1 /FULL/PATH/TO/CAV001_R1.fq.gz \
	-2 /FULL/PATH/TO/CAV001_R2.fq.gz \
	-o /FULL/PATH/TO/OUTPUT_DIR/CAV001

date
