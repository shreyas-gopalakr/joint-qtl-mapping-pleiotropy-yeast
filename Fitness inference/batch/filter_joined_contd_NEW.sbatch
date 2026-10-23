#! /bin/bash -e
#SBATCH -p shared # Partition to submit to (comma separated)
#SBATCH -J filter_contd_NEW # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-20:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o tmp/filter_contd_NEW_%A.out
#SBATCH -e tmp/filter_contd_NEW_%A.err

. ./vars.sh

bc1=($(awk 'BEGIN {FS = "\t"} {print $1}' data/complete_addresses_5thresh_bayes.txt))

bc2=($(awk 'BEGIN {FS = "\t"} {print $2}' data/complete_addresses_5thresh_bayes.txt))

psd_bc=$(tail -n 1 ${scratch}/filtered_counts_${seq_pool}.txt | awk 'BEGIN {FS = "\t"} {print $1}')

for i in "${!bc1[@]}" 
do if [[ ${bc1[i]}'_'${bc2[i]} == ${psd_bc} ]]
then psd_idx=$i
fi 
done

let arr_size=${#bc1[@]}-1
let start_idx=${psd_idx}+1
echo ${seq_pool}
for i in $(seq ${start_idx} ${arr_size})
do if [[ $(grep -P ${bc1[i]}'_'${bc2[i]} ${scratch}/joined_counts_${seq_pool}.txt) ]]; 
	then grep -P ${bc1[i]}'_'${bc2[i]} ${scratch}/joined_counts_${seq_pool}.txt >> ${scratch}/filtered_counts_${seq_pool}.txt #only grepping reads with barcode pairs found in BBQ segs
	echo $i
fi
done
