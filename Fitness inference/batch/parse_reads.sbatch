#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J parse # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-06:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 8000 # Memory in MB
#SBATCH -o tmp/parse_%A.out
#SBATCH -e tmp/parse_%A.err

module load centos6/0.0.1-fasrc01
module load python/3.4.1-fasrc01
source activate my_env

. ./vars.sh

#scratch=/n/holyscratch01/desai_lab/sugoshreyas/2021_8_26_sepsis_sequencing_pilot
#outdir=/n/desai_lab/users/sugoshreyas/SEPSIS_sequencing_pilot/secondary_data

# Clean reads from 4 lanes and join them 
bash clean_reads.sh data/fastq/SEPSIS/${sample_name}_S*_L001_R1_001.fastq.gz data/fastq/SEPSIS/${sample_name}_S*_L001_R2_001.fastq.gz > ${scratch}/cleaned_${sample_name}.txt
bash clean_reads.sh data/fastq/SEPSIS/${sample_name}_S*_L002_R1_001.fastq.gz data/fastq/SEPSIS/${sample_name}_S*_L002_R2_001.fastq.gz >> ${scratch}/cleaned_${sample_name}.txt
bash clean_reads.sh data/fastq/SEPSIS/${sample_name}_S*_L003_R1_001.fastq.gz data/fastq/SEPSIS/${sample_name}_S*_L003_R2_001.fastq.gz >> ${scratch}/cleaned_${sample_name}.txt
bash clean_reads.sh data/fastq/SEPSIS/${sample_name}_S*_L004_R1_001.fastq.gz data/fastq/SEPSIS/${sample_name}_S*_L004_R2_001.fastq.gz >> ${scratch}/cleaned_${sample_name}.txt

# Parse reads. Uses regular expression to split read into BC1, BC2, row and col indices (as numbers). P5 and P7 indices also here as numbers.
# How does parse know to wait until clean is done? Is clean just super fast?
cat ${scratch}/cleaned_${sample_name}.txt | python read_parser.py > ${scratch}/parsed_${sample_name}.txt
