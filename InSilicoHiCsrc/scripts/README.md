# In-silico Hi-C pipeline

This directory contains the code implementing the in-silico Hi-C workflow used in:

“Understanding the physical processes behind DNA–DNA proximity ligation assays”.

Each subdirectory in this folder corresponds to one stage of the pipeline.

The pipeline is organized into five stages:

## 1. FIRE minimization
Input chromatin structures are preprocessed and energy-minimized using the FIRE algorithm.

## 2. MD simulation with crosslinking and digestion
For each minimized structure, molecular dynamics simulations are performed using Langevin dynamics. Crosslinking to a protein matrix and enzymatic digestion of DNA are applied prior to time evolution. This stage produces MD trajectories and per-trajectory metadata (`.npz` files) describing crosslinked beads and fragment ends.

## 3. Ligation analysis
MD trajectories and per-trajectory metadata are analyzed to identify ligation events. A ligation event is recorded when two fragment ends come within a specified distance threshold and have not previously ligated, and when a probabilistic ligation criterion is satisfied. This stage produces one ligation log (`.npz` file) per trajectory.

## 4. Aggregation
Per-trajectory ligation logs are aggregated to generate final outputs, including ligation maps, contact maps, and ligation/contact probabilities as functions of genomic and Euclidean distance.

## 5. Plotting
Final aggregated `.npz` outputs are used to generate the figures included in the manuscript.
