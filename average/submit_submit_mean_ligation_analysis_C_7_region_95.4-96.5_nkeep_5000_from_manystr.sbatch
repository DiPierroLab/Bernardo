#!/bin/bash  
#SBATCH --job-name=submit_mean_ligation_analysis_C_7_region_95.4-96.5_nkeep_5000_from_manystr 
#SBATCH --output=submit_mean_ligation_analysis_C_7_region_95.4-96.5_nkeep_5000_from_manystr.out 
#SBATCH --error=submit_mean_ligation_analysis_C_7_region_95.4-96.5_nkeep_5000_from_manystr.err 
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --exclusive
#SBATCH --partition=ctbp

source ~/.bashrc
conda activate hic
python mean_first_time_ligation_analysis_C_7_region_95.4-96.5_nkeep_5000_from_manystr.py
#conda deactivate
