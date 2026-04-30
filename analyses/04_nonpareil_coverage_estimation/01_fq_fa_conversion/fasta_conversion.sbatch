#!/bin/bash
#SBATCH --job-name=fasta_conversion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -A <proj>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>



python3 ./fasta_conversion.py
