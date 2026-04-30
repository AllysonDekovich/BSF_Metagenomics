#!/bin/bash
#SBATCH --job-name=prokka_DM_T7_bac
#SBATCH --array=1-11
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>
#SBATCH --mem=16G

# set up job array
infile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" T0_bac_samples.txt)

prefix=$(basename "$infile" | sed 's/.fa*//g')

prokka \
--outdir "./output/${prefix}" \
--prefix "${prefix}" \
--kingdom Bacteria \
--gcode 11 \
--metagenome \
--mincontiglen 1000 \
--cpus 8 \
"${infile}"
