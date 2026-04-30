#!/bin/bash
#SBATCH --job-name=DM_T3_metabat2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


metabat2 \
-i DM_T3_final.contigs.fa \
-o ./DM_T3_bins \
-a DM_T3_coverage_metabat2.tsv \
-t 16
