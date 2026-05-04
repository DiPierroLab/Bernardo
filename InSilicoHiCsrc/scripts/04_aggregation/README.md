# Aggregation

This folder contains scripts for the aggregation stage of the in-silico Hi-C pipeline.

In this stage, per-trajectory outputs are combined into ensemble-averaged observables that form the final results used in the manuscript.

---

## Overview

Starting from outputs of the ligation-analysis stage (Stage 03) and MD trajectories (Stage 02), this step:

* aggregates ligation events across many independent trajectories
* constructs ensemble-averaged ligation maps
* computes ligation probabilities as functions of genomic and Euclidean distance
* computes trajectory-derived contact maps directly from MD configurations
* accumulates statistics of distances, displacements, and time evolution
* computes native-ensemble reference observables from minimized structures

The result is a set of aggregated datasets that can be directly compared with experimental Hi-C data.

---

## File overview

### Ligation aggregation script

`aggregate_ligation_statistics.py`

* Aggregates per-trajectory ligation-analysis outputs (`.npz` files)
* Computes:

  * ensemble-averaged ligation maps
  * ligation probability vs genomic distance
  * ligation probability vs Euclidean distance
  * cumulative ligation statistics and histories

---

### Time-resolved contact-map script

`compute_time_resolved_contact_maps.py`

* Computes average contact maps directly from MD trajectories at selected time frames
* Also computes:

  * average distance maps
  * contact probability vs genomic distance
  * mean-square displacement (MSD) over time
* Provides a trajectory-based **ground-truth reference** independent of the ligation process

---

### Native ensemble aggregation script

`compute_native_ensemble_maps.py`

* Computes average distance maps, contact maps, and contact probability curves from an ensemble of minimized native chromatin structures
* Reads minimized `.gsd` configurations and computes pairwise Euclidean distances
* Defines contacts using a distance threshold (`rc`)
* Averages over a subset of structures (`nsamples`)
* Produces the **native-ensemble baseline** used for comparison with crosslinked and ligated structures

Outputs include:

* average distance map
* average contact map
* contact probability vs genomic distance
* number of structures used (`nf`)

---

### Submission scripts

* `submit_aggregate_ligation_statistics.sbatch`
* `submit_time_resolved_contact_maps.sbatch`
* `submit_native_ensemble_maps.sbatch`

Example Slurm submission scripts for running the aggregation tasks on an HPC system.  
Users should adapt partition, environment, and resource settings.

---

## Inputs

This stage uses:

* Ligation-analysis outputs from Stage 03:

  * `new_firsttime_ligation_map_*.npz`

* MD trajectory data from Stage 02:

  * `chromosome_crosslinked_fixed_chopped_*.gsd`

* Metadata files:

  * `chopping_and_ligation_data_*.npz`

* Minimized native structures:

  * `minimized_*.gsd`

---

## Outputs

The aggregation step produces compressed `.npz` files containing:

### Ligation-derived observables

* ensemble-averaged ligation maps
* ligation probability vs genomic distance
* ligation probability vs Euclidean distance
* cumulative ligation histories

### Trajectory-derived observables

* average contact maps at different time points
* average distance maps
* contact probability vs genomic distance
* mean-square displacement (MSD)

### Native-ensemble observables

* average contact maps
* average distance maps
* contact probability vs genomic distance

These outputs are used to generate the figures in the manuscript and to compare simulated ligation signals with the underlying chromatin organization.

---

## Notes on reproducibility

* Aggregation operates over large ensembles of trajectories and structures
* Only final aggregated outputs are included in the repository
* Intermediate per-trajectory files can be regenerated using earlier stages
* Paths and cluster-specific settings have been generalized for portability
* Native-ensemble observables are computed from a subset of available minimized structures (`nsamples`)

---

## Role in the pipeline

This stage converts raw simulation outputs into the final statistical observables:

* **Ligation-based observables** emulate experimental Hi-C readouts
* **Trajectory-based observables** provide a ground-truth reference
* **Native-ensemble observables** define the baseline chromatin organization prior to crosslinking

Together, they enable a direct comparison between simulated proximity ligation data and the underlying chromatin structure.