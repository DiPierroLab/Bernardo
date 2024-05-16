#!/bin/bash  
#SBATCH --job-name=degraded_distance_contact_map 
#SBATCH --output=degraded_distance_contact_map.out 
#SBATCH --error=degraded_distance_contact_map.err 
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --exclusive
#SBATCH --partition=ctbp

source ~/.bashrc
conda activate hic
python degraded_distance_contact_map.py

#conda deactivate
