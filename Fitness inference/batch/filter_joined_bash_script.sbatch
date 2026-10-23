#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J filter # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-01:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o tmp/filter_%A.out
#SBATCH -e tmp/filter_%A.err

. ./vars.sh

bc1=($(awk 'BEGIN {FS = "\t"} {print $1}' data/complete_addresses_5thresh_bayes.txt))

bc2=($(awk 'BEGIN {FS = "\t"} {print $2}' data/complete_addresses_5thresh_bayes.txt))

for i in "${!bc1[@]}" 
do if [[ $(grep -P ${bc1[i]}'_'${bc2[i]} ${scratch}/joined_counts_${seq_pool}.txt) ]]; 
	then grep -P ${bc1[i]}'_'${bc2[i]} ${scratch}/joined_counts_${seq_pool}.txt >> ${scratch}/filtered_counts_${seq_pool}.txt #only grepping reads with barcode pairs found in BBQ segs
fi
done
