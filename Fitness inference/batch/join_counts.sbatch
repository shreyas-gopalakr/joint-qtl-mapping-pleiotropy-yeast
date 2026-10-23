#! /bin/bash -e
#SBATCH -p shared # Partition to submit to (comma separated)
#SBATCH -J reform # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-05:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/reformat_%A.out
#SBATCH -e tmp/reformat_%A.err

. ./vars.sh

module load centos6/0.0.1-fasrc01
module load python/3.4.1-fasrc01
source activate my_env

python join_counts.py ${seq_pool} # generates ${scratch}/joined_counts_${seq_pool}.txt
