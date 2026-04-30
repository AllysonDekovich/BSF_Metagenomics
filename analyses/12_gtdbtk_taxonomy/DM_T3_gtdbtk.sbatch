#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


gtdbtk classify_wf \
--genome_dir ./DM_T3 \
--extension fa \
--skip_ani_screen \
--out_dir ./DM_T3/output \
--cpus 8
