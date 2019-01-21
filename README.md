# lodestone

Lodestone implementation on Rivanna 

This pipeline orchestrates bioinformatics operations on a single isolate sequence data

Workflow Steps:
	QC --> Taxonomic assignment --> Variant calling --> De-novo assembly: spades

***

## Dependencies
	- java/1.8.0
	- python/2.7
	- perl/5.24.0
	- fastqc/0.11.5
	- trimgalore/0.4.5
	- cutadapt/1.16
	- trimmomatic/0.36
	- kaiju/1.6.3
	- mash/2.1
	- edirect/1.0.0
	- bwa/0.7.17
	- samtools/1.19
	- sambamba/0.6.8
	- picard/2.18.5
	- htslib/1.9
	- bcftools/1.9
	- vcftools/0.1.15
	- vcflib/8.22
	- spades/3.11.1
	- quast/5.0.1
	- bedtools/2.26.0

***

## Usage

### lodestone.sh

This script performs QC, followed by reference identification using kaiju/mash. The reference genome is downloaded from NCBI, followed by variant calling. The pipeline will also perform de-novo assembly using spades. 
 
```
[~]$ sh lodestone.sh [options]

OPTIONS:
        -h      Show this message
        -c      CAVID
        -1      FULL/PATH/TO/ForwardRawReads (*_R1.fq.gz)
        -2      FULL/PATH/TO/ReverseRawReads (*_R2.fq.gz)
        -o      FULL/PATH/TO/OutDir

## Example command 
[~]$ sh lodestone.sh -c CAV001 -1 CAV001_R1.fq.gz -2 CAV002_R2.fq.gz -o CAV001_lodestone

```

### lodestone_w_ref.sh

Use this script if you wish to call variants against a known reference genome.   
This script performs QC, followed by variant calling. The pipeline will also perform de-novo assembly using spades. 

```
[~]$ sh lodestone_w_ref.sh [options]

OPTIONS:
        -h      Show this message
        -c      CAVID
        -1      FULL/PATH/TO/ForwardRawReads (*_R1.fq.gz)
        -2      FULL/PATH/TO/ReverseRawReads (*_R2.fq.gz)
        -r      FULL/PATH/TO/ref.fa (Reference genome FASTA)
        -o      FULL/PATH/TO/OutDir

## Example command 
[~]$ sh lodestone_w_ref.sh -c CAV001 -1 CAV001_R1.fq.gz -2 CAV002_R2.fq.gz -r ref.fa -o CAV001_lodestone

```

### SLURM submission scripts

Example SLURM job submission scripts are included
	- `submit-lodestone.slurm.sh` : job submission for single isolate data
	- `submit-lodestone.slurm-array.sh`	: batch submission for multiple isolates using slurm array




