#!/bin/bash
#SBATCH --job-name=compute_native_ensemble_maps
#SBATCH --output=compute_native_ensemble_maps.out
#SBATCH --error=compute_native_ensemble_maps.err
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --partition=partition_name

# Load user shell configuration, if needed
source ~/.bashrc

# Activate conda environment
conda activate hic

# Run script for computing average contact maps, distance maps,
# and contact probabilities from the native ensemble
python compute_native_ensemble_maps.py