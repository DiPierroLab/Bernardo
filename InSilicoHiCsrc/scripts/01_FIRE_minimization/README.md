# FIRE minimization

This folder contains scripts for the initial energy minimization stage of the in-silico Hi-C pipeline.

In this stage, chromatin structures are relaxed using the FIRE (Fast Inertial Relaxation Engine) algorithm prior to molecular dynamics simulations.

---

## Overview

Starting from input chromatin configurations (e.g., PDB trajectories), this stage:

* loads individual structures
* constructs a polymer model representation
* performs FIRE energy minimization
* outputs minimized configurations as `.gsd` files

These minimized structures are used as input for the subsequent MD crosslinking/digestion stage.

---

## File overview

### Template script

`minimize_chromosome_configuration_source.py`

* Template Python script used to generate all minimization runs
* Contains placeholder variables of the form `@variable`
* Implements:

  * structure loading and centering
  * bond and angle assignment
  * FIRE energy minimization
  * short Langevin relaxation
  * output of minimized structure (`.gsd`)

---

### Generator script

`run_many_minimize.sh`

* Generates one Python script per structure
* Replaces placeholders in the template script
* Creates job lists and Slurm submission scripts
* Used to launch large-scale minimization runs

Paths and cluster-specific settings have been generalized for portability.

---

### Slurm submission template

`submit_chr_all_source`

* Template used to generate Slurm batch submission scripts
* Uses GNU parallel to execute jobs from a job list
* Placeholders (`@jobname`, `@joblist`, `@ntasks`) are filled by the generator script

---

### Example generated script

`minimize_energy_*.py`

* Example of a fully instantiated minimization script
* Corresponds to one specific structure and random seed
* Included for illustration only; paths should be adapted before execution

---

### Example submission script

`submit_chr_*.sbatch`

* Example Slurm submission script generated from the template
* Illustrates how jobs are grouped and submitted

---

## Notes on reproducibility

* Production runs generated many per-structure scripts and job files
* These are not included, as they can be regenerated from the template and generator scripts
* Paths and cluster-specific settings have been generalized
* Execution of this stage requires an HPC environment with HOOMD-blue and Slurm

---

## Outputs

Each minimization produces:

* `.gsd` files containing minimized chromatin configurations

These outputs are used as input for the MD crosslinking/digestion stage.

