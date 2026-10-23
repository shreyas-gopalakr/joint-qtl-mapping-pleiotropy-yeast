#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J number_of_correctly_indexed_reads # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/number_of_correctly_indexed_reads_%A.out
#SBATCH -e tmp/number_of_correctly_indexed_reads_%A.err

sample_names=($(awk 'BEGIN {FS = ","} {print $3}' sample_sheet_sepsis_novaseq_run1.csv)) #has the wrong kcl env names, but that's fine for our purposes here
sample_name_list=($(awk 'BEGIN {FS = ","} {print $15}' data/primers_novaseq_run1_final.csv))
p5_id_list=($(awk 'BEGIN {FS = ","} {print $13}' data/primers_novaseq_run1_final.csv))
p7_id_list=($(awk 'BEGIN {FS = ","} {print $14}' data/primers_novaseq_run1_final.csv))
col_id_list=($(awk 'BEGIN {FS = ","} {print $11}' data/primers_novaseq_run1_final.csv))
row_id_list=($(awk 'BEGIN {FS = ","} {print $12}' data/primers_novaseq_run1_final.csv))
for i in {0..1044}
do
	echo \> ${sample_names[$i]}
	for j in {1..1045}
	do if [[ ${sample_name_list[j]} == ${sample_names[$i]} ]]
	then idx=$j
	fi
	done
	p5_id=${p5_id_list[idx]}
	p7_id=${p7_id_list[idx]}
	col_id=${col_id_list[idx]}
	row_id=${row_id_list[idx]}
	correct_index_reads=$(grep -P ${p5_id}'\t'${p7_id}'\t'${col_id}'\t'${row_id}'\t' scratch/corrected_${sample_names[$i]}.txt | wc -l)
	total_index_reads=$(cat scratch/corrected_${sample_names[$i]}.txt | wc -l)
	echo -e ${sample_names[$i]}'\t'${correct_index_reads}'\t'${total_index_reads} >> number_of_correctly_indexed_reads.txt
done

