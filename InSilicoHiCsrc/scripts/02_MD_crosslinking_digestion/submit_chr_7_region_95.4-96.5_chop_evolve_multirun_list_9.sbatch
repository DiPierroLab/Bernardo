#!/bin/bash  
#
# Example Slurm submission script generated from submit_chr_all_source
# by run_chop_evolve_from_many_str.sh.
#
# This file illustrates how MD crosslinking/digestion jobs were bundled
# and submitted in production runs.
#
# It is included for documentation only; partition, environment, and
# resource settings should be adapted before execution.

#SBATCH --job-name=chr_7_region_95.4-96.5_chop_evolve_multirun_list_9 
#SBATCH --output=chr_7_region_95.4-96.5_chop_evolve_multirun_list_9.out 
#SBATCH --error=chr_7_region_95.4-96.5_chop_evolve_multirun_list_9.err 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=115  
#SBATCH --time=100:00:00
#SBATCH --exclusive
#SBATCH --partition=partition_name

# Load user shell configuration (if needed)
source ~/.bashrc

# Activate conda environment
conda activate 

# Run job list using GNU parallel
cat chr_7_region_95.4-96.5_chop_evolve_multirun_list_9.txt | path/to/parallel

