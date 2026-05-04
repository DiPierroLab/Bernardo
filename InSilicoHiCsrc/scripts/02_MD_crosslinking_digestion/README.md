# MD crosslinking and digestion

This folder contains scripts for the molecular dynamics (MD) stage of the in-silico Hi-C pipeline, where chromatin structures undergo crosslinking, enzymatic digestion (bond removal), and time evolution.

Each subdirectory corresponds to one stage of the pipeline; this folder implements the MD simulation stage.

---

## File overview

### 1. Template script

```text
chromosome_fix_chop_evolve_from_singlestr_source.py
```

This is the **template script** used to generate all production MD simulation scripts.

* Contains placeholder variables of the form `@variable`
* Implements:

  * bond removal (DNA digestion)
  * identification of fragment ends
  * assignment of evolving vs fixed particles (crosslinking)
  * Langevin dynamics evolution
  * output of trajectory (`.gsd`) and metadata (`.npz`)

---

### 2. Generator script

```text
run_chop_evolve_from_many_str.sh
```

This shell script:

* iterates over input chromatin structures
* assigns random seeds and parameters
* replaces placeholders in the template script
* generates **one Python script per structure and per run**
* creates job lists and Slurm submission scripts

This script was used to launch large-scale simulation campaigns.

---

### 3. Slurm submission template

```text
submit_chr_all_source
```

Template used to generate Slurm batch submission scripts.

* placeholders (`@jobname`, `@joblist`, `@ntasks`) are filled by the generator script
* executes jobs using GNU parallel

---

### 4. Example generated script

```text
run_chr_chop_evolve_*.py
```

Example of a **fully instantiated simulation script** produced from the template.

* corresponds to one specific structure and one random seed
* included for illustration only
* uses generic placeholder paths for portability

---

### 5. Example submission script

```text
submit_chr_*.sbatch
```

Example Slurm submission script generated from the template.

* shows how jobs are bundled and submitted to the scheduler

---

## Notes on reproducibility

* Production runs generated **thousands of per-structure scripts and job files**
* These are **not included** because they can be regenerated from the template and generator scripts
* Paths and cluster-specific settings have been generalized for portability
* Execution of this stage requires an HPC environment with HOOMD-blue and Slurm

---

## Outputs

Each simulation produces:

* `.gsd` trajectory files (MD evolution)
* `.npz` metadata files containing:

  * ligating fragment ends
  * removed bonds and angles
  * fixed particle indices

These outputs are used as input for downstream ligation analysis

