#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J index # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-05:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/index_%A.out
#SBATCH -e tmp/index_%A.err

. ./vars.sh

sample_name_list=($(awk 'BEGIN {FS = ","} {print $15}' data/primers_novaseq_run1_final.csv)) #this file needs to have sample name and the 4 primers
for i in {1..1045}
do if [[ ${sample_name_list[i]} == ${sample_name} ]]
then idx=$i
fi
done

p5_id_list=($(awk 'BEGIN {FS = ","} {print $13}' data/primers_novaseq_run1_final.csv))
p7_id_list=($(awk 'BEGIN {FS = ","} {print $14}' data/primers_novaseq_run1_final.csv))
col_id_list=($(awk 'BEGIN {FS = ","} {print $11}' data/primers_novaseq_run1_final.csv))
row_id_list=($(awk 'BEGIN {FS = ","} {print $12}' data/primers_novaseq_run1_final.csv))

p5_id=${p5_id_list[idx]}
p7_id=${p7_id_list[idx]}
col_id=${col_id_list[idx]}
row_id=${row_id_list[idx]}

grep -P ${p5_id}'\t'${p7_id}'\t'${col_id}'\t'${row_id}'\t' ${scratch}/corrected_${sample_name}.txt | sort | uniq > ${scratch}/indexed_${sample_name}.txt #only grepping reads with the correct multiplexing indices
