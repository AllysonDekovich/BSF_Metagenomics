#!/bin/bash
#SBATCH --job-name=nonpareil_redundancy_kmer
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu
#SBATCH --array=1-9


# initiate job array
# filenames.txt contains all files that need to be submitted to array
infile=$(sed -n -e "${SLURM_ARRAY_TASK_ID} p" filenames.txt)
#echo $infile

outfile=$(basename "$infile" | sed 's/.fastq.gz//')
#echo $outfile

nonpareil -s "$infile" -T kmer -k 15 -f fastq -b "$outfile"
