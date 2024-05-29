# README file for DNA degradation simulations

Included in this folder are source code files for generating Python scripts that simulate the degradation of an ensemble of native structures representative of lymphoblastoid cells over a 1.1Mb region of chromosome 7 (95.4 to 96.5 Mbp) at a nucleosome resolution (200bp), comprised of 5500 beads each. 

## Ensemble of original structures

For details on the ensemble of structures, see:
Harris, H.L., Gu, H., Olshansky, M., Wang, A., Farabella, I., Eliaz, Y., Kalluchi, A., Krishna, A., Jacobs, M., Cauer, G. and Pham, M., 2023. Chromatin alternates between A and B compartments at kilobase scale for subgenic organization. Nature communications, 14(1), p.3303.
https://www.nature.com/articles/s41467-023-38429-1

The ensemble of structures was generated using NuChroM.
Source code for NuChroM can be accessed at: 
https://github.com/DiPierroLab/NuChroM

## Relaxation of original structures

The native structures from the ensemble of original structures are first subjected to a relaxation scheme that uses the FIRE algorithm to minimize the energy of the samples.
The Hamiltonian includes the following forces: FENE bonds, harmonic angle potentials and truncated Lennard-Jones potentials.
After the FIRE algorithm has converged, the structures are further evolved in time with Langevin dynamics for the purpose of thermalization.
The final output of this initial relaxation code for an individual structure is a '.gsd' file that contains the energy-optimized structure which is ready to be used in the DNA degradation simulations.

The shell script 'run_many_minimize.sh' takes in the file 'minimize_chromosome_configuration_source.py' as a template to generate a set of individual Python scripts, one for each structure in the ensemble.
The shell script also submits these Python scripts as batch jobs to a computing cluster under SLURM, using the file 'submit_chr_all_source' as a template for generating and submitting individual batch scripts.

## Simulation of DNA degradation

The DNA degradation process takes the individual, energy-optimized, native structures from the .gsd files, evolves them in time with Langevin dynamics and chops a total of number of bonds 'ncut' at random periodically, where the period between successive cuts is given by a total number of time steps 'nsteps'. 
This periodic process of bond cutting is repeated for a total number of degradation cycles equal to 'ncycles'.

The interaction in the Langevin dynamics includes: FENE bonds, harmonic angle potentials and truncated Lennard-Jones potentials.

The particular protocol used in the simulations is such that all bonds are chopped at the very beginning of the simulation and the only interaction that remains is that of excluded volume (i.e. Lennard-Jones).
That is to say: all 5499 bonds are cut at t=0. This specific protocol is achieved by setting ncut=5499 and ncycles=1. In this way, the chromosomes bonds are fully degraded at the very start and only once cycle of cutting is required, thereby producing a gas of 5500 nucleosomes that is allowed to evolve in time.

The shell script 'run_degrade_from_many_gas.sh' takes in the file 'chromosome_degrade_from_singlestr_source_gas.py' as a template to generate a set of individual Python scripts, one for each structure in the ensemble. 
The shell script also submits these Python scripts as batch jobs to a computing cluster under SLURM, using the file 'submit_chr_all_source' as a template for generating and submitting individual batch scripts.

The Python scripts require the use of the molecular dynamics package HOOMD version 2.9.7.
The final outputs of a individual Python script includes '.gsd' files that store the trajectories of the nucleosomes of the degraded structure as well as an '.npz' file with data on which bonds were cut and at what time.

## Calculation of average distance map, average contact map, mean-square displacements and contact probability.

With the ensemble of molecular trajectories stored as '.gsd' files and the auxiliary '.npz' files, the calculation of average distance maps, average contact maps, contact probabilities and average mean-square displacements of the nucleosomes follows by running the Python script 'degraded_distance_contact_map.py'.

In this Python script, one can choose the individual frames of the trajectories on which one wants to calculate average distance and contact maps and contact probabilities.
The individual frames chosen be specified in the list 'iframes'.
'nsamples=5000' is the total number of trajectories to sample over, 'size=5500' is the number of nucleosomes in a sample, 'rc=1.5' specifies the cut-off distance (15 nm) below which two nucleosomes are considered to be in contact. 
Since only one degradation cycle is considered, 'ncycle=1' and 'icycle=0'.

The final outputs with the average distance map, average contact map, mean-square displacements and contact probability are then stored in an '.npz' file, which can then be opened and read with a Jupyter notebook like the .ipynb file included for further data analysis and for visualization.

## Python package requirements

The current version of the code runs on an Anaconda environment. 
The list of packages installed in this conda environment can be found in the file 'requirements_conda.txt'
