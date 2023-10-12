#!/bin/bash 
#SBATCH -p        # put comma separated list of available partitions here 
#SBATCH -n 1                    # Number of cores
#SBATCH -t 0-3:00               # Runtime in days-hours:minutes 
#SBATCH --mem 6000              # Memory in MB
#SBATCH -J rmunfix               # job name 
#SBATCH -o rmunfix.%A.out        # File to which standard out will be written 
#SBATCH -e rmunvix.%A.err        # File to which standard err will be written 
##SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
##SBATCH --mail-user=<PUT YOUR EMAIL ADDRESS HERE>  # Email to which notifications will be sent

r1=$1 # rCorrector corrected left (R1)fastq file
r2=$2 # rCorrector corrected right (R2) fastq file
sample_name=$3 # sample name used for outfiles prefix
echo "input read pair files are: $r1 and $r2""

FilterUncorrectabledPEfastq.py -1 $r1 -2 $r2 -s $sample $sample_name

