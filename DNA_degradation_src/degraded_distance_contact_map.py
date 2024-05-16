#!/usr/bin/env python
# coding: utf-8

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
new_output_dir = '/scratch/b.zubillagaherrera/woolly_init_chop_2/degraded_distance_contact_map/5500_degcuts_RMSD_200_v2'
input_dir   = '/scratch/b.zubillagaherrera/woolly_init_chop_2/degrade_ncut_5500_long'

nstart = 1
nfiles = 20
nsims = 1
gamma = 1
size = 5500
rc = 1.5
nsamples = 5000 #Total number of samples to average over
Nframes = 5000 #Total number of time frames
icycle = 0 #Degradation used from whose last frame the Hi-C protocol was initiated
ncycle = 1
ncut = 5499 #Number of cuts per cycle during degradation
nsteps = 5000000

#iframes = [0, 99, 249, 499, 749, 999, 1249, 2499, 4999] #Must include frame 0 as first element of the list for MSD calculations
iframes = [0, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213] #Must include frame 0 as first element of the list for MSD calculations

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
            fn_base = 'C_'+str(c)+'_region_'+str(r)+'_ncut_'+str(ncut)+'_nsteps_'+str(nsteps)+'_ncycle_'+str(ncycle)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(isim)
            
            degradation_data_file = input_dir+'/degradation_data_'+fn_base+'.npz'
            degraded_str_file     = input_dir+'/degraded_'+fn_base+'.gsd'
            
            if (os.path.exists(degradation_data_file)) and (os.path.exists(degraded_str_file)) and (np.min(n_files)<nsamples):

                #reading input file containing output of trajectory with degraded chromosome over long time
                f = gsd.fl.open(degraded_str_file,'rb')

                #reading file that contains details of bonds and angles kept after degradation, as well as particles at free ends
                g = np.load(degradation_data_file,allow_pickle=True)

                fin_bonds_record=g['fin_bonds_record']
                ncycles=g['ncycles']
                ncut=g['ncut']

                #Pick the individual frame of the trajectory corresponding to the individual cycle chosen
                #to use as initial structure/snapshort for this calculation.
                #Frame picked is the very last frame of the chosen cycle for maximal diffusion.
                Nframes_in_file = gsd.hoomd.HOOMDTrajectory(f).__len__() #total number of frames in gsd trajectory file
                NframesPerCycle = int(Nframes/ncycles)  #number of frames per cycle
                iframe = NframesPerCycle*(icycle+1) - 1 #number of the frame chosen as initial structure for this calculation.

                print('frame number ', iframe, ' cycle number ', icycle, ' no. frames per cycle ', NframesPerCycle, ' no. frames ', Nframes, ' no. cycles', ncycles)
                print('desired frame ', iframe, 'no. of frames in .gsd file', Nframes_in_file)
                print('desired cycle ', icycle, 'no. of cycles in .npz file', len(fin_bonds_record))
                
                nf += 1
                print('n ', nf, degraded_str_file, ' found')     
                
                isnap = 0
                
                for iframe in iframes:
                    
                    if (iframe <= Nframes_in_file-1) and (n_files[isnap]<nsamples):
                    
                        print(iframe)
                    
                        #creating hoomd snapshot with final frame at particular degradation cycle in consideration
                        config = hoomd.data.gsd_snapshot(degraded_str_file,frame=iframe)
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

new_fn_base = 'C_'+str(c)+'_region_'+str(r)+'_ncut_'+str(ncut)+'_nsteps_'+str(nsteps)+'_ncycle_'+str(ncycle)+'_icycle_'+str(icycle)+'_rc_'+str(rc)

'C_'+str(c)+'_region_'+str(r)+'_ncut_'+str(ncut)+'_nsteps_'+str(nsteps)+'_ncycle_'+str(ncycle)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(isim)
output_file = new_output_dir+'/degraded_distance_contact_map_'+new_fn_base+'_diff.npz'

np.savez_compressed(output_file,iframes=iframes,nsnapshots=nsnapshots,distance_map_average=distance_map_average,contact_map_average=contact_map_average,contact_prob_average=contact_prob_average, nf=nf, n_files=n_files, mean_square_disp=mean_square_disp)

print(str(nf)+' files are used')




