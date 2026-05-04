#!/usr/bin/env python
# coding: utf-8

# Compute time-resolved contact maps from MD trajectories.
#
# This script averages contact maps, distance maps, contact probabilities,
# and mean-square displacements over an ensemble of MD trajectories at
# selected time frames.
#
# These trajectory-derived contact maps provide a reference for comparison
# against simulated ligation maps.
#
# Inputs:
#   - MD trajectories from stage 02 (.gsd)
#   - Per-trajectory metadata files (chopping_and_ligation_data_*.npz)
#
# Outputs:
#   - A compressed .npz file containing average contact maps, distance maps,
#     contact probabilities, sampled frame indices, and MSD values.
#
# Paths and parameters should be adapted before execution.

import hoomd, hoomd.md
import mdtraj as md
import mbuild as mb
import numpy  as np
import random
#from joblib import Parallel, delayed
import multiprocessing
import os.path
import sys
import scipy as sp
from scipy.spatial import distance
#import matplotlib.pyplot as plt
#import matplotlib as mpl
import gsd
import gsd.fl, gsd.hoomd
import warnings
import copy

warnings.filterwarnings('ignore')

c = 7  # chromosome identity
r = '95.4-96.5'  # region identity
new_output_dir = '/path/to/aggregation_outputs'
input_dir = '/path/to/md_outputs'

nstart = 1
nfiles = 20
nsims = 1
gamma = 1
size = 5500
rc = 1.5
nsamples = 5000 #Total number of samples to average over
nkeep = 4500

ncut = 1000

iframes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 250, 500, 750, 1000, 1250, 2500, 4999] #Must include frame 0 as first element of the list for MSD calculations
nsnapshots = len(iframes)

contact_map_average  = np.zeros((nsnapshots, size, size))
distance_map_average = np.zeros((nsnapshots, size, size))

mean_square_disp = np.zeros(nsnapshots)

hoomd.context.initialize("");

nf = 0

n_files = np.zeros(nsnapshots)

for filenumber in range(1,nfiles+1):
    for str_no in range(0,999,1):    
        for isim in range(1,nsims+1):
            fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nkeep)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(isim)
            
            data_file = input_dir+'/chopping_and_ligation_data_'+fn_base+'.npz'
            str_file     = input_dir+'/chromosome_crosslinked_fixed_chopped_'+fn_base+'.gsd'
            
            if (os.path.exists(data_file)) and (os.path.exists(str_file)) and (np.min(n_files)<nsamples):

                #reading input file containing output of trajectory with degraded chromosome over long time
                f = gsd.fl.open(str_file,'rb')

                #Total number of frames in file
                Nframes_in_file = gsd.hoomd.HOOMDTrajectory(f).__len__() #total number of frames in gsd trajectory file


                print(f"Processing sample {nf+1}: {fn_base}")

                nf += 1
                   
                
                isnap = 0
                
                for iframe in iframes:
                    
                    if (iframe <= Nframes_in_file-1) and (n_files[isnap]<nsamples):
                    
                        print('frame no. ', iframe)
                    
                        #creating hoomd snapshot with final frame at particular degradation cycle in consideration
                        config = hoomd.data.gsd_snapshot(str_file,frame=iframe)
                        pos = config.particles.position
                        distance_map = distance.cdist(pos, pos, 'euclidean')
                    
                        size = config.particles.N
                        contact_map  = np.zeros((size, size))
                        distance_map_thresholded = rc - distance_map
                        contact_map  = (np.sign(distance_map_thresholded)+1.)/2.
    
                        distance_map_average[isnap] += distance_map
                        contact_map_average[isnap]  += contact_map
            
                        #Mean squared displacement - initial frame is taken as reference for displacement
                        if (iframe == 0):
                            initial_frame_pos = copy.deepcopy(pos)
                        disp_vec = pos - initial_frame_pos
                        msq = np.square(disp_vec)
                        msd = np.average(msq)
                        
                        mean_square_disp[isnap] += msd 
                        
                        n_files[isnap] += 1
            
                    isnap += 1
                
                
print(np.max(contact_map_average), n_files)
    
#Calculation of ligation map and ligation probability vs genomic distance without cap
dmax = size - 1
init = 1
for isnap in range(nsnapshots):
    if n_files[isnap]>0:
        contact_map_average[isnap] /= n_files[isnap]

contact_prob_average = np.zeros((nsnapshots, dmax))
for isnap in range(nsnapshots):
    for i in range(dmax):
        contact_prob_average[isnap][i] =  np.mean(np.diagonal(contact_map_average[isnap], offset=(i+init)))    
    
#Calculation of final average distance map
for isnap in range(nsnapshots):
    if n_files[isnap]>0:
        distance_map_average[isnap] /= n_files[isnap]
        
#Mean squared displacement for each time frame averaged over all samples
for isnap in range(nsnapshots):
    if n_files[isnap]>0:
        mean_square_disp[isnap] /= n_files[isnap]

new_fn_base = 'C_'+str(c)+'_region_'+str(r)+'_ncut_'+str(ncut)

output_file = new_output_dir+'/evolved_distance_contact_map_'+new_fn_base+'.npz'

np.savez_compressed(output_file,iframes=iframes,nsnapshots=nsnapshots,distance_map_average=distance_map_average,contact_map_average=contact_map_average,contact_prob_average=contact_prob_average, nf=nf, n_files=n_files, mean_square_disp=mean_square_disp)

print(str(nf)+' files are used')




