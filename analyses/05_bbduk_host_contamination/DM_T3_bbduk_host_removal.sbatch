#!/bin/bash
#SBATCH --job-name=bbduk_DM_T3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --array=1-3

# Only select lines with _R1_001.fastq.gz and pick the correct one by SLURM_ARRAY_TASK_ID
# DM_T3_filenames.txt contains all samples that will be submitted to the array
infile1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" DM_T3_filenames.txt)
echo "R1: $infile1"

# Replace _R1_001.fastq.gz with _R2_001.fastq.gz
infile2=$(basename "$infile1" | sed 's/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz/')
echo "R2: $infile2"

#create outfile basename only
outfile=$(basename "$infile1" | sed 's/_R1_trimmed.fastq.gz//')
echo "outfile: $outfile"

# run bbduk

java -ea -Xmx80g -Xms80g -cp /nfs/home/adekovic/mambaforge/envs/bbmap/opt/bbmap-39.33-0/current/ jgi.BBDuk \
in1=$infile1 \
in2=$infile2 \
out1=${outfile}_R1_host_trimmed.fastq.gz \
out2=${outfile}_R2_host_trimmed.fastq.gz \
ref=GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz \ #downloaded from NCBI (see manuscript for link)
k=25 \
hdist=0 \
ktrim=r \
outm1=bos_taurus_R1.fastq.gz \
outm2=bos_taurus_R2.fastq.gz \
tpe \
stats={$outfile}_alignment_stats.txt
