#!/bin/bash
#SBATCH --job-name=coverm_T0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


# run coverm on all samples, will average them out downstream via R
coverm contig \
--coupled DM_T7_r1_Owings_S13_L001_R1_host_trimmed.fastq.gz DM_T7_r1_Owings_S13_L001_R2_host_trimmed.fastq.gz \
DM_T7_r2_Owings_S14_L001_R1_host_trimmed.fastq.gz DM_T7_r2_Owings_S14_L001_R2_host_trimmed.fastq.gz \
DM_T7_r3_Owings_S15_L001_R1_host_trimmed.fastq.gz DM_T7_r3_Owings_S15_L001_R2_host_trimmed.fastq.gz \
--reference DM_T7_final.contigs.fa \
--methods tpm \
--min-read-percent-identity 95 \
--min-read-aligned-percent 75 \
--output-file DM_T7_tpm.tsv
