#!/bin/bash
#SBATCH --job-name=DM_C_T0_maxbin2_binning
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>


run_MaxBin.pl -contig DM_C_T0_final.contigs.fa -abund DM_C_T0_coverage_maxbin2.tsv \
-min_contig_length 1500 -thread 16 -out DM_C_T0_bins
