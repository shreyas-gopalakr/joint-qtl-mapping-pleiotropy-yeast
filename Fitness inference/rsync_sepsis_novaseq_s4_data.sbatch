#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J rsync_sepsis_novaseq_s4_data # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o tmp/rsync_sepsis_novaseq_s4_data_%A.out
#SBATCH -e tmp/rsync_sepsis_novaseq_s4_data_%A.err

rsync -ah --progress /n/ngsdata/220919_A00794_0697_AH5G5VDSX5_i5revcomp/ /n/desai_lab/users/sugoshreyas/SEPSIS/novaseq_run1/raw_data
