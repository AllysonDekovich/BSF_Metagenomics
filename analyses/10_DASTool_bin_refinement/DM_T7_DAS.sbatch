#!/bin/bash
#SBATCH --job-name=DM_T7_DAS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


DAS_Tool \
-i ./DM_T7_maxbin2_contig2bin.tsv,./DM_T7_metabat2_contig2bin.tsv \
-l maxbin2,metabat \
-c ./DM_T7_final.contigs.fa \
-o ./DM_T7_refined_bins \
-t 30
