#!/bin/bash
#SBATCH --job-name=DM_T7_checkM
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adekovic@vols.utk.edu


checkm2 predict \
--threads 30 \
--input ./DM_T7_bins \
--output-directory ./DM_T7_checked_bins \
-x fa \
--database_path ./CheckM2_database/uniref100.KO.1.dmnd
