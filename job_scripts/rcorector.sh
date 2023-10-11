#!/bin/bash 
#SBATCH -p        # put comma-separated list of available partisions here 
#SBATCH -n 16                   # Number of cores
#SBATCH -t 23:00:00               # Runtime in hours:minutes:secods
#SBATCH --mem 48000               # Memory in MB
#SBATCH -J rcorrector               # job name 
#SBATCH -o rcorrector.%A.out        # File to which standard out will be written 
#SBATCH -e rcorrector.%A.err        # File to which standard err will be written 
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=<PUT YOUR EMAIL ADDRESS HERE>  # Email to which notifications will be sent

"""
For this script to initialize a conda environment, a version of python that supports
anaconda or mamba will need to be in PATH. In a case where the HPC environment has python
available as a loadable module, such as the Harvard Cannon cluster, 
it will simply require adding: module load python.
"""

source activate rcorrector

"""
R1 and R2 are the first and second command line arguments that follow "sbatch rcorrector.sh
R1 is a comma-separated list of the left (R1) fastq files, and R2 is a similar list for the right (R2) reads
"""

R1=$1
R2=$2

run_rcorrector.pl -t 16 -1 $R1 -2 $R2
