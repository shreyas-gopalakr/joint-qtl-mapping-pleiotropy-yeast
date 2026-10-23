#! /bin/bash -e
#SBATCH -p shared # Partition to submit to (comma separated)
#SBATCH -J infer_rep1 # Job name
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:30 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o tmp/infer_rep1_%A.out
#SBATCH -e tmp/infer_rep1_%A.err

module load gcc/8.2.0-fasrc01
module load R/3.3.3-fasrc01

. ./vars.sh

#setting priors to 0 means that the script ignores these priors. The Rscript calculates the means and variances from each round to pass on to the next round.

prior_mean=0
prior_var=0
echo ${env_label}
for i in {1..3}
do
	echo it $i mean = ${prior_mean} var = ${prior_var}
	/n/home00/nnguyenba/scripts/PCG_BFA_with_prior_5 ${prior_mean} ${prior_var} ${scratch}/uniq_filtered_counts_header_added/header_added_${env_label}.txt > ${scratch}/fitness_${env_label}_it${i}.txt
	prior_mean=$(Rscript calculate_fitness_mean.R ${scratch}/fitness_${env_label}_it${i}.txt)
	prior_var=$(Rscript calculate_fitness_var.R ${scratch}/fitness_${env_label}_it${i}.txt)
done

rsync ${scratch}/fitness_${env_label}_it${i}.txt ${outdir}/fitness/rep1/fitness_${env_label}_it${i}.txt
