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
warnings.filterwarnings('ignore')

nc = @nkeep
c = @chr  # chromosome identity
r = '@region'  # region identity
output_dir = '@output_dir'
new_output_dir = '@new_output_dir'
nstart = @nstart
nfiles = @nfiles
nsims = @nsims
gamma = @gamma
size = @size
rc = @rc
nsamples = @nsamples #Total number of samples to average over

lig_probs = np.array([1.0,0.1,0.01,0.001,0.0001])
n_prob = len(lig_probs)
ligation_map_global = np.zeros((n_prob,size, size))
dmax = size - 1
init = 1
ligation_prob_global = np.zeros((n_prob,dmax))
nf = 0
tot_nlig = np.zeros(n_prob)
tot_slig = np.zeros(n_prob)
n_lseg = 0
n_pairs = 0
dl = 5
dc = 5

init_contact_map_final = np.zeros((size,size))
initial_distance_final = np.zeros((size,size))
firstfile = 0

rthresh   = rc*np.array([10, 7.5, 5, 2.5, 1])
n_rthresh = len(rthresh)
init_contact_map_final_2 = np.zeros((len(rthresh),size,size))

dx = 0.2
x_edges = np.arange(dx,10+dx,dx)
all_fij_global = np.zeros((n_prob,len(x_edges)-1))
fij_global = np.zeros((n_prob,len(x_edges)-1))
p_ligation1 = np.zeros((n_prob,len(x_edges)-1))
p_ligation2 = np.zeros((n_prob,len(x_edges)-1))
all_fij = np.zeros((n_prob,len(x_edges)-1))
all_Dij = np.zeros((n_prob,len(x_edges)-1))
#length_ligation_2dhist = []
#length_ligation_events = []
max_length = 0

for filenumber in range(1,nfiles+1):
    for str_no in range(0,999,1):
        for isim in range(1,nsims+1):
            fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(isim)
            ligationmap_file = new_output_dir+'/new_firsttime_ligation_map_'+fn_base+'_rc_'+str(rc)+'.npz'
            
            if (os.path.exists(ligationmap_file)) and (nf<nsamples):
                
                
                data = np.load(ligationmap_file,allow_pickle=True)
                initial_distance_mat = data["initial_distance_mat"]
                initial_distance_final += initial_distance_mat
                init_contact_map_global = data["init_contact_map_global"]
                init_contact_map_final += init_contact_map_global
                
                #Calculation of original contact maps from original distance maps
                init_contact_map_global_2 = np.zeros((len(rthresh),size, size))
                for irt in range(n_rthresh):
                    distance_mat_thresholded        = rthresh[irt] - initial_distance_mat
                    init_contact_map_global_2[irt]  = (np.sign(distance_mat_thresholded)+1.)/2.
                init_contact_map_final_2 += init_contact_map_global_2

                
                g = np.load(output_dir+'/chopping_and_ligation_data_'+fn_base+'.npz')
                
                #print(fn_base)
                
                ligating_ends = g["ligating_ends"]
                ligating_ends = np.sort(np.append(ligating_ends,[0,size-1],0))
                n_sp = len(ligating_ends)
                n_spr = int(n_sp*(n_sp-1)/2)

                #print('lig ends',ligating_ends)
            
                data = np.load(ligationmap_file,allow_pickle=True) 
                
                Pcycle_L   = data["Pcycle_L"]
                cycle_bins = data["cycle_bins"]
                Nlong_L    = data["Nlong_L"]
                if firstfile == 0:
                    Nlong_L_final  = np.copy(Nlong_L)
                    Pcycle_L_final = np.copy(Pcycle_L)
                    firstfile == 1
                else:
                    Pcycle_L_final += Pcycle_L
                    Nlong_L_final  += Nlong_L
                
                ligation_map = data["ligation_map"]
                self_ligation_events = data["self_ligation_events"]
                long_segments = data["long_segments"]
                fij_global = data["fij_global"]
                Dij = data["Dij"]
                for mm in range(n_prob):
                    fij = fij_global[mm]*Dij
                    all_fij[mm] += fij
                    all_Dij[mm] += Dij
                    all_fij_global[mm] += fij_global[mm]
                    #bb = [ len(x) for x in long_segments ]
                    #cc = sp.sum(self_ligation_events[mm].tolist(),axis=0)
                    #cc = cc.tolist()[0]
                    #l_max = np.max(bb)
                    #if l_max > max_length:
                    #    max_length = l_max
                    #c_max = np.max(cc)
                    #l_edges = np.arange(0,l_max+dl,dl)
                    #c_edges = np.arange(0,c_max+dc,dc)
                    #H,xedges,yedges = np.histogram2d(bb,cc,bins=(l_edges,c_edges))
                    #length_ligation_2dhist.append(H)
                    #ee = np.sign(cc)
                    #length_ligation_events.extend(np.where(ee > 0)[0])
                    tot_nlig[mm] += np.sum(ligation_map[mm])
                    #tot_slig[mm] += np.sum(np.sign(cc))
                    ligation_map_global[mm] += ligation_map[mm]
                    p_ligation1[mm] += fij_global[mm]
                n_lseg += len(long_segments)
                n_pairs += n_spr
                nf += 1
                
                print('sample number '+str(nf)+' done')
                
    print('file number '+str(filenumber)+' done')

init = 1
for mm in range(n_prob):
    ligation_map_global[mm] /= nf
    for i in range(dmax):
        ligation_prob_global[mm][i] =  np.mean(np.diagonal(ligation_map_global[mm], offset=(i+init)))

#ll_edges = np.arange(0,max_length+dl,dl)
#h_length, ll_edges = np.histogram(length_ligation_events,bins=ll_edges)

    p_ligation1[mm] = all_fij_global[mm]/nf
    p_ligation2[mm] = all_fij[mm]/all_Dij[mm]
    
    #print('p_ligation1',p_ligation1)
    #print('p_ligation2',p_ligation2)

init_contact_map_final /= nf
initial_distance_final /= nf

init_contact_map_final_2 /=nf

#Cyclization probabilities
Ncycle_L_final = np.copy(Pcycle_L_final)
for i in range(n_prob):
    Pcycle_L_final[i] = np.divide(Pcycle_L_final[i], Nlong_L_final) 
    
ligation_distances = x_edges

new_fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_rc_'+str(rc)
mean_ligationmap_file = new_output_dir+'/combined_ft_lig_data_'+new_fn_base+'.npz'
#mean_ligationmap_file = output_dir+'/combined_first_time_ligation_map_'+new_fn_base+'.npz'

np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,init_contact_map_final_2=init_contact_map_final_2,rthresh=rthresh,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_L_final=Nlong_L_final,cycle_bins=cycle_bins)

#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_final=Nlong_final,cycle_bins=cycle_bins)



#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances)
#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=[tot_slig,n_lseg],length_ligation_2dhist=length_ligation_2dhist,dl=dl,dc=dc,nf=nf,h_length=h_length,ll_edges=ll_edges,ligation_distance_prob=ligation_distance_prob,ligation_distances=ligation_distances)
print(str(nf)+' files are used')
print(str(tot_nlig)+' ligation events out of '+str(n_pairs)+' open end pairs')
#print(str(tot_slig)+' total self ligation events out of '+str(n_lseg)+' long segments')


