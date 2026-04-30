#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --array=1-9

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
# filenames.txt is a list of all the files to be submitted (*.fastq.gz)
infile1=$(grep "_R1_001.fastq.gz" filenames.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
#echo "R1: $infile1"

# Replace _R1_001.fastq.gz with _R2_001.fastq.gz
infile2=$(basename "$infile1" | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
#echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_001.fastq.gz//')
#echo "outfile: $outfile"

# run fastp
fastp \
--in1 "$infile1" \
--in2 "$infile2" \
--out1 ./output/"${outfile}_R1_trimmed.fastq.gz" \
--out2 ./output/"${outfile}_R2_trimmed.fastq.gz" \
--adapter_sequence TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
--adapter_sequence_r2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
--detect_adapter_for_pe \
--qualified_quality_phred 25 \
--html ./output/"${outfile}".html \
--thread 4
