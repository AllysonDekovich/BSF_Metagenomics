#!/bin/bash
#SBATCH --job-name=fastqc_pretrimmed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<user email>

for file in *fastq.gz
do
	fastqc -t 6 -o QC_output_pretrim $file
done
