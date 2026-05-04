#!/usr/bin/env python
# coding: utf-8

#  compute_native_ensemble_maps.py 
#  Computes average distance maps, contact maps, and contact probability curves from an ensemble of minimized native chromatin structures. 
#  The script reads minimized `.gsd` structures, computes pairwise Euclidean distances, thresholds contacts using the distance 'rc', 
#  averages over `nsamples` structures, and saves the results to a compressed `.npz` file.
#
#  Output keys:
#  - `distance_map_average`
#  - `contact_map_average`
#  - `contact_prob_average`
#  - `nf`
#
#  This script generates the native-ensemble baseline used for comparison with crosslinked and ligated structures.

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
warnings.filterwarnings('ignore')

c = 7  # chromosome identity
r = '95.4-96.5'  # region identity
new_output_dir = '/path/to/native_ensemble'
output_dir   = '/path/to/minimized_structures'

nstart = 1
nfiles = 20
nsims = 1
gamma = 1
size = 5500
rc = 1.5
nsamples = 5 #Total number of samples to average over

contact_map_average  = np.zeros((size, size))
distance_map_average = np.zeros((size, size))

hoomd.context.initialize("");

nf = 0

for filenumber in range(1,nfiles+1):
    for str_no in range(0,999,1):
        for isim in range(1,nsims+1):
            fn_base = 'C_'+str(c)+'_region_'+str(r)+'_file_'+str(filenumber)+'_str_'+str(str_no)
            minimized_file = output_dir+'/minimized_'+fn_base+'_minimized.gsd'
                        
            if (os.path.exists(minimized_file)) and (nf<nsamples):     

                #reading input file containing trajectory
                f = gsd.fl.open(minimized_file,'rb')
                
                config = hoomd.data.gsd_snapshot(minimized_file,frame=f.nframes-1)
                pos = config.particles.position
                distance_map = distance.cdist(pos, pos, 'euclidean')
    
                size = config.particles.N
                contact_map  = np.zeros((size, size))
                distance_map_thresholded = rc - distance_map
                contact_map  = (np.sign(distance_map_thresholded)+1.)/2.

                distance_map_average += distance_map
                contact_map_average  += contact_map
                    
                nf += 1
                print('n ', nf, minimized_file, ' found')

    
#Calculation of ligation map and ligation probability vs genomic distance without cap
dmax = size - 1
init = 1
contact_map_average /= nf
contact_prob_average = np.zeros(dmax)
for i in range(dmax):
    contact_prob_average[i] =  np.mean(np.diagonal(contact_map_average, offset=(i+init)))    
    
#Calculation of final average distance map

distance_map_average /= nf

new_fn_base = 'C_'+str(c)+'_region_'+str(r)

output_file = new_output_dir+'/distance_contact_map_'+new_fn_base+'.npz'

np.savez_compressed(output_file,distance_map_average=distance_map_average,contact_map_average=contact_map_average,contact_prob_average=contact_prob_average, nf=nf)

print(str(nf)+' files are used')




