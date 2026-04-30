#!/bin/bash
#SBATCH --job-name=megahit_T3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>



megahit \
-1 DM_T3_r1_Owings_S16_L001_R1_host_trimmed.fastq.gz,DM_T3_r3_Owings_S17_L001_R1_host_trimmed.fastq.gz,DM_T3_r4_Owings_S18_L001_R1_host_trimmed.fastq.gz \
-2 DM_T3_r1_Owings_S16_L001_R2_host_trimmed.fastq.gz,DM_T3_r3_Owings_S17_L001_R2_host_trimmed.fastq.gz,DM_T3_r4_Owings_S18_L001_R2_host_trimmed.fastq.gz \
-t 12 \
-m 260000000000 \
--k-list 21,29,39,49,59,69,79,89,99,109,119,129,139,149,159,169,179,189,199 \
--min-contig-len 1000 \
-o DM_T3
