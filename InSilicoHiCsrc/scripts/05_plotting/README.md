# Plotting scripts

This folder contains Jupyter notebooks used to generate figures and analyze results from the in-silico Hi-C simulations.

These notebooks operate on processed datasets (trimmed `.npz` files) and reproduce the main trends and selected supplementary figures presented in the manuscript.

## Notebooks

- `plot_crosslinking_time_evolution.ipynb`  
  Explores how chromatin structure evolves over time for different crosslinking strengths.  
  Includes contact probability curves, mean-square displacement, and representative contact maps.

- `plot_ligation_statistics.ipynb`  
  Analyzes global ligation statistics and distributions derived from simulated Hi-C experiments.

- `plot_ligation_crosslinking.ipynb`  
  Examines how ligation outcomes depend on crosslinking strength.

- `plot_ligation_single_cell_statistics.ipynb`  
  Computes and visualizes variability and statistics at the single-cell level.

- `plot_KR_normalization_and_PCA.ipynb`  
  Applies KR normalization and principal component analysis to contact maps.

## Data requirements

All notebooks expect input data in: `../../data/processed/`

In particular:
- `crosslinking_time_evolution/` for time-resolved observables  
- additional processed `.npz` files depending on the analysis  

The `.npz` datasets required for these notebooks are **not included in this repository** due to size constraints.

They are available on Zenodo:  
https://doi.org/10.5281/zenodo.19682913

All notebooks in this folder are designed to operate on those datasets.

## Notes

- The provided datasets are **trimmed** to reduce size while preserving all key observables.  
- Full grids of contact and distance maps are not included; only representative cases are provided where needed.  
- These notebooks are intended to **illustrate and reproduce key results**, not to regenerate every intermediate dataset from scratch.  

## Usage

Run the notebooks in order of interest. No strict execution order is required, as each notebook is self-contained and loads its own data.

