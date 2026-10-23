#! /bin/bash -e
#SBATCH -p shared # Partition to submit to (comma separated)
#SBATCH -J copy # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o tmp/copy_%A.out
#SBATCH -e tmp/copy_%A.err

mv /n/desai_lab/users/sugoshreyas/SEPSIS/novaseq_run1/scratch/filtered_counts_* /n/desai_lab/users/sugoshreyas/SEPSIS/novaseq_run1/scratch/filtered_bcs_incomplete  
