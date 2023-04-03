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

import math

nc = 5000
c = 7  # chromosome identity
r = '95.4-96.5'  # region identity
new_output_dir = '/scratch/b.zubillagaherrera/time_test/launch_analysis/nnsamples_5000'
output_dir   = '/scratch/b.zubillagaherrera/time_test/rc1.5'
output_dir_2 = '/scratch/b.zubillagaherrera/time_test/fix_chop_evolve'

nstart = 1
nfiles = 20
nsims = 1
gamma = 1
size = 5500
rc = 1.5
nsamples = 3 #5000 #Maximum total number of files/samples to average over

ncounts = 1000 #Total number of Hi-C map counts to consider

lig_probs = np.array([1.0,0.1,0.01,0.001])
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


n_hits  = np.zeros(n_prob)
n_files = np.zeros(n_prob)
count_thresh = 0

#I changed the order to sample more evenly and fairly across files
#Counts will be accumulated until a threshold (ncounts) is reached, irrespective of the number of files used to average over.
#This is to guarantee the same number of counts in each ligation map
for isim in range(1,nsims+1):
    for str_no in range(0,999,1):
        for filenumber in range(1,nfiles+1):
            fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(isim)
            ligationmap_file = output_dir+'/new_firsttime_ligation_map_'+fn_base+'_rc_'+str(rc)+'.npz'
            
            #print(ligationmap_file)
            
            if (os.path.exists(ligationmap_file) and (nf<nsamples)):
                
                
                data = np.load(ligationmap_file,allow_pickle=True)
                initial_distance_mat = data["initial_distance_mat"]
                initial_distance_final += initial_distance_mat
                init_contact_map_global = data["init_contact_map_global"]
                init_contact_map_final += init_contact_map_global
                
                print(ligationmap_file)
                
                #Calculation of original contact maps from original distance maps
                init_contact_map_global_2 = np.zeros((len(rthresh),size, size))
                for irt in range(n_rthresh):
                    distance_mat_thresholded        = rthresh[irt] - initial_distance_mat
                    init_contact_map_global_2[irt]  = (np.sign(distance_mat_thresholded)+1.)/2.
                init_contact_map_final_2 += init_contact_map_global_2

                
                g = np.load(output_dir_2+'/chopping_and_ligation_data_'+fn_base+'.npz')
                
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
                    firstfile = 1
                else:
                    Pcycle_L_final += Pcycle_L
                    Nlong_L_final  += Nlong_L
                
                ligation_map = data["ligation_map"]
                self_ligation_events = data["self_ligation_events"]
                long_segments = data["long_segments"]
                fij_global = data["fij_global"]
                Dij = data["Dij"]
                
                for mm in range(n_prob):
                    if (n_hits[mm] < ncounts):
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
                        
                        n_hits[mm]  += np.sum(ligation_map[mm])
                        n_files[mm] += 1
                        count_thresh = min(n_hits)
                    
                n_lseg += len(long_segments)
                n_pairs += n_spr
                nf += 1
                
                print('sample number '+str(nf)+' done')
                
                print(n_hits, count_thresh, nf, n_files, tot_nlig)
                
    print('file number '+str(filenumber)+' done')

print(n_hits, count_thresh, n_files, tot_nlig)

print(all_fij[0],all_Dij[0],all_fij_global[0])

print(np.sum(all_fij[0][3:]))

#sys.exit()


#We now trim the ligation vs. Euclidean distance counts to guarantee that thet have exactly the same
#desired total number of ligations, which is equal to ncounts. We do this by removing hits from the maps at distribution at random,
#according to the weights of the dist. itself, until they all have the same final number of counts ncounts.
    
n_fij_init = np.zeros(n_prob)

for mm in range(n_prob):
    
    n_bins = len(all_fij[mm])
    
    #Exclude bins that contain nan
    bin_start = 0
    for i in range(n_bins):
        if (math.isnan(all_fij[mm][i])==True):
            bin_start +=1
    
    #Number of ligation events to get rid off
    n_fij_init[mm] = np.sum(all_fij[mm][bin_start:])
    n_decimate = int(n_fij_init[mm]) - int(ncounts)
    print(np.sum(all_fij[mm][bin_start:]), ncounts, type(ncounts), n_decimate, type(ncounts))
    
    #List with bin numbers
    bin_numbers = np.arange(bin_start,n_bins,1)
    print(bin_numbers)
    print(n_bins)

    #We sample the lig. prob. vs. Eucledian distance distribution and pick the ligations to get rid off.
    pick = random.choices(bin_numbers, weights = all_fij[mm][bin_start:], k=1000*n_decimate)

    npicks = 0
    ite = 0
    while(npicks < n_decimate):
        xx = int(pick[ite])
        occupancy = int(all_fij[mm][xx])
        
        print('xx', xx)
        print('iteration',npicks)
        
        if (occupancy>0):
            print(ite, pick[ite], xx, occupancy)
            print('pre', all_fij[mm][xx])
            all_fij[mm][xx] -= 1
            npicks += 1
            print('post', all_fij[mm][xx])
            
        ite +=1
        
    
    
    print('p = ',lig_probs[mm],' lig prob vs distance trimmed from ',tot_nlig[mm],' to ',np.sum(all_fij[mm][bin_start:]),' ligations')     
    
    sys.exit()

#Now we proceed to trim off the excesses in order to guarantee that all the ligation maps have exactly the same
#desired total number of ligations, which is equal to ncounts. We do this by removing hits from the maps at random,
#according to the weights of each matrix element, until they all have the same final number of counts ncounts.

#First we create a 1d array with the coordinates of each matrix element.
coord_list = []
for i_locus in range(size):
    for j_locus in range(size):
        coord_list.append([i_locus, j_locus])
        #print(coord_list)       
        
for mm in range(n_prob):
    
    #Now we flatten the ligation maps into 1d arrays whose coordinates correspond to those in coord_list.
    flat_lig_map = ligation_map_global[mm].flatten()
    
    #Number of ligation events to get rid off
    n_decimate = int(n_hits[mm]) - int(ncounts)
    print(n_hits[mm], type(n_hits[mm]), ncounts, type(ncounts), n_decimate, type(ncounts))
    
    print(len(flat_lig_map))
    print(len(coord_list))
    
    ##We sample the ligation map according to the corresponding weights and pick the ligations to get rid off.
    pick = random.choices(coord_list, weights = flat_lig_map, k=1000*n_decimate)
    
    npicks = 0
    ite = 0
    while(npicks < n_decimate):
        xx = int(pick[ite][0])
        yy = int(pick[ite][1])
        occupancy = int(ligation_map_global[mm][xx][yy])
        
        #print('iteration',npicks)
        
        if (occupancy>0):
            #print(ite, pick[ite], xx, yy, occupancy)
            #print('pre', ligation_map_global[mm][xx][yy])
            ligation_map_global[mm][xx][yy] -= 1
            npicks += 1
            #print('post', ligation_map_global[mm][xx][yy])
            
        ite +=1
        
    n_hits[mm] = np.sum(ligation_map_global[mm])
    
    print('p = ',lig_probs[mm],' trimmed from ',tot_nlig[mm],' to ',n_hits[mm],' ligations') 
    

init = 1
for mm in range(n_prob):
    ligation_map_global[mm] /= n_files[mm] #nf
    for i in range(dmax):
        ligation_prob_global[mm][i] =  np.mean(np.diagonal(ligation_map_global[mm], offset=(i+init)))

#ll_edges = np.arange(0,max_length+dl,dl)
#h_length, ll_edges = np.histogram(length_ligation_events,bins=ll_edges)

    p_ligation1[mm] = all_fij_global[mm]/n_files[mm] #nf
    p_ligation2[mm] = all_fij[mm]/all_Dij[mm]
    
    #print('p_ligation1',p_ligation1)
    #print('p_ligation2',p_ligation2)

init_contact_map_final /= nf
initial_distance_final /= nf

init_contact_map_final_2 /= nf

#Cyclization probabilities
Ncycle_L_final = np.copy(Pcycle_L_final)
for i in range(n_prob):
    Pcycle_L_final[i] = np.divide(Pcycle_L_final[i], Nlong_L_final) 
    
ligation_distances = x_edges

new_fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_rc_'+str(rc)
mean_ligationmap_file = new_output_dir+'/combined_ft_lig_data_same_counts'+new_fn_base+'.npz'
#mean_ligationmap_file = output_dir+'/combined_first_time_ligation_map_'+new_fn_base+'.npz'

np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances,init_contact_map_final=init_contact_map_final,init_contact_map_final_2=init_contact_map_final_2,rthresh=rthresh,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_L_final=Nlong_L_final,cycle_bins=cycle_bins,n_hits=n_hits)

#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,init_contact_map_final_2=init_contact_map_final_2,rthresh=rthresh,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_L_final=Nlong_L_final,cycle_bins=cycle_bins)

#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_final=Nlong_final,cycle_bins=cycle_bins)



#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances)
#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=[tot_slig,n_lseg],length_ligation_2dhist=length_ligation_2dhist,dl=dl,dc=dc,nf=nf,h_length=h_length,ll_edges=ll_edges,ligation_distance_prob=ligation_distance_prob,ligation_distances=ligation_distances)
#print(str(nf)+' files are used')
print(str(n_files)+'files are used')
print(str(tot_nlig)+' ligation events out of '+str(n_pairs)+' open end pairs')
#print(str(tot_slig)+' total self ligation events out of '+str(n_lseg)+' long segments')


