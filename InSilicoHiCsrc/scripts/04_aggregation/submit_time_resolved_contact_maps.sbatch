#!/bin/bash  
#
# Example Slurm submission script for time-resolved contact-map aggregation.
#
# This script runs the trajectory-based contact-map averaging step, which
# computes average contact maps, distance maps, contact probabilities, and
# mean-square displacements from MD trajectories.
#
# It is included for documentation purposes; users should adapt partition,
# environment, walltime, and paths to their own HPC system.

#SBATCH --job-name=time_resolved_contact_maps
#SBATCH --output=time_resolved_contact_maps.out
#SBATCH --error=time_resolved_contact_maps.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=partition_name

source ~/.bashrc
conda activate hic

python compute_time_resolved_contact_maps.py
