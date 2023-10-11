#!/bin/bash 
#SBATCH -p        # put comma separated list of available partitions here 
#SBATCH -n 1                    # Number of cores
#SBATCH -t 0-3:00               # Runtime in days-hours:minutes 
#SBATCH --mem 6000              # Memory in MB
#SBATCH -J FastQC               # job name 
#SBATCH -o FastQC.%A.out        # File to which standard out will be written 
#SBATCH -e FastQC.%A.err        # File to which standard err will be written 
##SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
##SBATCH --mail-user=<PUT YOUR EMAIL ADDRESS HERE>  # Email to which notifications will be sent

"""
For this script to initialize a conda environment, a version of python that supports
anaconda or mamba will need to be in PATH. In a case where the HPC environment has python
available as a loadable module, such as the Harvard Cannon cluster, 
it will simply require adding: module load python.
"""

source activate fastqc
infile=$1 # $1 represent the first (and in this case only) command line argument supplied to the script
echo "infile is $infile"

fastqc --outdir `pwd`/fastqc $infile

