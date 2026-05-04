#!/bin/bash  
#
# Example Slurm submission script generated from submit_chr_all_source
# by run_many_minimize.sh.
#
# This file illustrates the structure of a generated submission script
# for a specific batch of FIRE minimization jobs. It is included for
# documentation only; paths, partition, environment, and job list names
# should be adapted before execution.

#SBATCH --job-name=chr_7_region_95.4-96.5_minimize_list_100 
#SBATCH --output=chr_7_region_95.4-96.5_minimize_list_100.out 
#SBATCH --error=chr_7_region_95.4-96.5_minimize_list_100.err 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20  
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --partition=partition_name

# Load user shell configuration (if needed)
source ~/.bashrc

# Activate conda environment
conda activate hic

# Run job list using GNU parallel
cat chr_7_region_95.4-96.5_minimize_list_100.txt | path/to/parallel
