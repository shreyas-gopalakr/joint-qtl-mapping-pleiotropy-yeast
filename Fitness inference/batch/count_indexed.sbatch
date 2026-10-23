#!/bin/bash
#
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J count # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-06:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/fa_counts_%A.out
#SBATCH -e tmp/fa_counts_%A.err

. ./vars.sh

cut -f 5,6 ${scratch}/indexed_${sample_name}.txt | sort | uniq -c | sort -nr | awk "{print \$2,\$3,\$1}" OFS="\\t" > ${scratch}/counts_${sample_name}.txt

# "\\t does the same job as "\t". However, when you do echo -e, it also prints a new line before the tab. Don't know why this is.
