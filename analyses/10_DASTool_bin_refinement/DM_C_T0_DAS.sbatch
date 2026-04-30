#!/bin/bash
#SBATCH --job-name=DM_C_T0_DAS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


DAS_Tool \
-i ./DM_C_T0_maxbin2_contig2bin.tsv,./DM_C_T0_metabat2_contig2bin.tsv \
-l maxbin2,metabat \
-c ./DM_C_T0_final.contigs.fa \
-o ./DM_C_T0_refined_bins \
-t 30
