#!/bin/bash
#SBATCH --job-name=fairy_coverage_DM_T7
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


# run fairy to calculate coverage on indexed files; create 2 coverage files per timepoint (for maxbin and metabat2)
fairy coverage DM_T7_final.contigs.fa ./DM_T7_indexed_files/*.bcsp -t 8 --maxbin-format -o DM_T7_coverage_maxbin2.tsv
fairy coverage DM_T7_final.contigs.fa ./DM_T7_indexed_files/*.bcsp -t 8 -o DM_T0_coverage_metabat2.tsv
