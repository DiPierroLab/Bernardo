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
import gsd
import gsd.fl, gsd.hoomd
import itertools
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings('ignore')

seed_num = @seed_num
random.seed(seed_num)
c = @chr  # chromosome identity
r = '@region'  # region identity
filenumber = @filenumber
str_no = @str
#run_no = @run
output_dir = '@output_dir'
new_output_dir = '@new_output_dir'
fn_base = '@fn_base'
rc = @rc
lig_probs = np.array([1.0,0.9,0.5,0.3,0.1,0.01])
n_prob = len(lig_probs)

input_dir = '@input_dir'
minimized_str_file = input_dir+'/minimized_C_'+str(c)+'_region_'+r+'_file_'+str(filenumber)+'_str_'+str(str_no)+'_minimized.gsd'
f0 = gsd.fl.open(minimized_str_file,'rb')
initial_config = gsd.hoomd.HOOMDTrajectory(f0).read_frame(f0.nframes-1)
pos = initial_config.particles.position
initial_distance_mat = distance.cdist(pos, pos, 'euclidean')

trajfile = output_dir+'/chromosome_crosslinked_fixed_chopped_'+fn_base+'.gsd'
contactmap_file = new_output_dir+'/new_firsttime_ligation_map_'+fn_base+'_rc_'+str(rc)+'.npz'

if not os.path.exists(contactmap_file):
    f = gsd.fl.open(trajfile,'rb')
    s0 = gsd.hoomd.HOOMDTrajectory(f).read_frame(0)
    nt = f.nframes
    size = s0.particles.N
    
    g = np.load(output_dir+'/chopping_and_ligation_data_'+fn_base+'.npz')
    ligating_ends = g["ligating_ends"]
    ligating_ends = np.sort(np.append(ligating_ends,[0,size-1],0))
    rbs = g["removed_bonds"]
    n_sp = len(ligating_ends)
    n_spr = int(n_sp*(n_sp-1)/2)
    ligating_pairs = []
    for pair in itertools.combinations(ligating_ends,2):
        ligating_pairs.append(pair)
    
    lp = np.array(ligating_pairs)
    e2 = 0
    segments = []
    long_segments = []
    for i in range(len(rbs)):
         e1 = rbs[i]
         seg1 = np.arange(e2,e1+1)
         e2 = e1+1
         segments.append(seg1)
         if len(seg1) >= 5:
             long_segments.append(seg1)
    
    n_seg = len(segments)
    n_lseg = len(long_segments)
    
    contact_map_global = np.zeros((size, size))
    #ligation_map_global = np.zeros((size, size))
    ligation_map_with_prob = np.zeros((n_prob,size,size))
    self_ligation_events = np.zeros((n_prob,nt,n_lseg))
    Ntotal = 0
    ligated_ends = []
    available_pairs = ligating_pairs
    ap = []
    n_avpr = np.zeros(n_prob)
    for i in range(n_prob):
        ap.append(np.array(available_pairs))
        n_avpr[i] = len(available_pairs)
    
    dx = 0.2
    x_edges = np.arange(dx,10+dx,dx)
    fij_global = np.zeros((n_prob,len(x_edges)-1))
    x_centers = 0.5*(x_edges[:-1]+x_edges[1:])
    for i in range(nt):
        s0 = gsd.hoomd.HOOMDTrajectory(f).read_frame(i)
        pos = s0.particles.position
        distance_mat = distance.cdist(pos, pos, 'euclidean')
        distance_mat_thresholded =  rc - distance_mat
        contact_map = (np.sign(distance_mat_thresholded)+1.)/2.
    
        contact_map_global += contact_map
        ligation_map = np.zeros((n_prob,size,size))
        for mm in range(n_prob):
            px = lig_probs[mm]
            for j in range(int(n_avpr[mm])):
                if contact_map[tuple(ap[mm][j])] > 0.:
                #print(j,available_pairs[j])
                    z = np.random.random()
                    if z <= px:
                        ligation_map[mm][tuple(ap[mm][j])] = 1

                        for m in range(n_lseg):
                            ancilla1 = [long_segments[m][0], long_segments[m][-1]]
                            ancilla2 = list(ap[mm][j])

                            ancilla1.sort()
                            ancilla2.sort()
                     
                            if ancilla1  == ancilla2:
                                self_ligation_events[mm][i][m] += 1
        
            ligation_map_with_prob[mm] += ligation_map[mm]
            nx,ny = np.where(ligation_map[mm] > 0.)
            nn = np.unique(np.append(nx,ny))
            for q in nn:
                ap[mm] = np.delete(ap[mm],np.where(ap[mm]==q)[0],axis=0)
            n_avpr[mm] = len(ap[mm]) #available_pairs)
            #print(i,n_avpr,ligation_map_global.max(),ligation_map.max())
    
        Ntotal += 1
        if i % 100 == 0:
            print("Reading frame {:} of {:}".format(i, f.nframes))
    
    ligatable_distance_mat = []
    for i in ligating_pairs:
        ligatable_distance_mat.append(initial_distance_mat[i])
    
    #sampling distance and contact prob
    self_lig_ev = []
    for i in range(n_prob):
        ligation_distance_mat = np.multiply(ligation_map_with_prob[i],initial_distance_mat)
        fij, x_edges = np.histogram(ligation_distance_mat.flat, bins=x_edges, density=False)
        Dij, x_edges = np.histogram(ligatable_distance_mat,bins=x_edges, density=False)
        fij_global[i] = fij/Dij
    
        self_lig_ev.append(csr_matrix(self_ligation_events[i]))
    #

    np.savez_compressed(contactmap_file,ligation_map=ligation_map_with_prob,Ntotal=Ntotal,self_ligation_events=self_lig_ev,long_segments=long_segments,Dij=Dij,fij_global=fij_global,x_edges=x_edges)
#
