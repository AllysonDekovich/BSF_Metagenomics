#!/bin/bash
#SBATCH --job-name=fairy_index_DM_C_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-3

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
infile1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" DM_C_T0_filenames.txt)
echo "R1: $infile1"

# create variable name for the second paired end read set
infile2=$(basename "$infile1" | sed 's/_R1_host_trimmed.fastq.gz/_R2_host_trimmed.fastq.gz/')
echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_host_trimmed.fastq.gz//')
echo "outfile: $outfile"

# run fairy to index short reads
fairy sketch -1 $infile1 -2 $infile2 -d DM_C_T0_indexed_files
