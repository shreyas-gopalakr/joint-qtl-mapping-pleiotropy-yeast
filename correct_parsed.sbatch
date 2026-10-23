#! /bin/bash -e
#SBATCH -p shared # Partition to submit to (comalta separated)
#SBATCH -J corr # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-04:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 50000 # Memory in MB
#SBATCH -o tmp/errcorr_%A.out
#SBATCH -e tmp/errcorr_%A.err

module load gcc/8.2.0-fasrc01

. ./vars.sh

/n/home00/nnguyenba/scripts/BBQ/error_correct_from_dict -bc 5 -err 3 -dict data/BY_dictionary.txt ${scratch}/parsed_${sample_name}.txt > ${scratch}/corrected_${sample_name}.tmp
/n/home00/nnguyenba/scripts/BBQ/error_correct_from_dict -bc 6 -err 3 -dict data/RM_dictionary.txt ${scratch}/corrected_${sample_name}.tmp > ${scratch}/corrected_${sample_name}.txt
rm ${scratch}/corrected_${sample_name}.tmp

#The BY and RM dictionary files are lists of BY barcodes and RM barcodes (30k BY and 15k RM barcodes)
#The error_correct_from_dict is a binary file.
