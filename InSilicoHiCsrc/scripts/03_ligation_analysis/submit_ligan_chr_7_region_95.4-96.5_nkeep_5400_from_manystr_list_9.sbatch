#!/bin/bash  
#
# Example Slurm submission script generated from submit_chr_all_source_v1
# by gen_firsttime_ligation_maps_from_manystr_v6.sh.
#
# This file illustrates how ligation-analysis jobs were bundled and submitted
# in production runs. It is included for documentation only; partition,
# environment, and paths should be adapted before execution.

#SBATCH --job-name=ligan_chr_7_region_95.4-96.5_nkeep_5400_from_manystr_list_9 
#SBATCH --output=ligan_chr_7_region_95.4-96.5_nkeep_5400_from_manystr_list_9.out 
#SBATCH --error=ligan_chr_7_region_95.4-96.5_nkeep_5400_from_manystr_list_9.err 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=115  
#SBATCH --time=168:00:00
#SBATCH --exclusive
#SBATCH --partition=partition_name


# Load user shell configuration (if needed)
source ~/.bashrc

# Activate conda environment
conda activate hic

# Run job list using GNU parallel
cat ligan_chr_7_region_95.4-96.5_nkeep_5400_from_manystr_list_9.txt | /path/to/parallel
