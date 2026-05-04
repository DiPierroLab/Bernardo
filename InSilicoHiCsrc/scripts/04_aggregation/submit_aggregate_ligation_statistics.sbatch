#!/bin/bash
#
# Example Slurm submission script for aggregation stage of the in-silico Hi-C pipeline.
#
# This script runs the aggregation of per-trajectory ligation outputs into
# ensemble-averaged observables.
#
# It is included for documentation purposes; users should adapt partition,
# environment, and paths to their own HPC system.

#SBATCH --job-name=aggregate_ligation_statistics
#SBATCH --output=aggregate_ligation_statistics.out
#SBATCH --error=aggregate_ligation_statistics.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=partition_name

# Load user shell configuration (if needed)
source ~/.bashrc

# Activate conda environment
conda activate hic

# Run aggregation script
python aggregate_ligation_statistics.py