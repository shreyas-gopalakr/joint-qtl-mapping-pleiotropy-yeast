#!/usr/bin/bash
#SBATCH --job-name=annotating_all_SNPs_cluster
#SBATCH --output=../shell_outputs/output_files/annotating_all_SNPs_cluster_%A_%a.out
#SBATCH --error=../shell_outputs/error_files/annotating_all_SNPs_cluster_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH -p dpetrov,normal
#SBATCH -c 1
#SBATCH --mem=8GB

module load matlab

#sbatch --array=0-41 annotating_all_SNPs_cluster.sbatch to run the code. Run it from current folder
matlab -nodisplay -nosplash -r "annotating_all_SNPs_cluster(${SLURM_ARRAY_TASK_ID}); exit"
