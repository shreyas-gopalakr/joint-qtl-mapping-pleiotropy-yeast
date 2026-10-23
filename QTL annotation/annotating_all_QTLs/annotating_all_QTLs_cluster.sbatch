#!/usr/bin/bash
#SBATCH --job-name=annotating_all_QTLs_cluster
#SBATCH --output=../shell_outputs/output_files/annotating_all_QTLs_cluster_%A.out
#SBATCH --error=../shell_outputs/error_files/annotating_all_QTLs_cluster_%A.err
#SBATCH --time=8:00:00
#SBATCH -p dpetrov,hns,normal
#SBATCH -c 1
#SBATCH --mem=8GB

module load matlab

#sbatch annotating_all_QTLs_cluster.sbatch to run the code. Run it from current folder
matlab -nodisplay -nosplash -r "annotating_all_QTLs; exit"
