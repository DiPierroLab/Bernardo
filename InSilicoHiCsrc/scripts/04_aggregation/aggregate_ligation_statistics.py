#!/usr/bin/env python
# coding: utf-8

# Aggregation script for in-silico Hi-C ligation analysis.
#
# This script aggregates per-trajectory ligation-analysis outputs (.npz files)
# into ensemble-averaged observables.
#
# Inputs:
#   - Per-trajectory ligation outputs from stage 03 (ligation analysis)
#   - Corresponding metadata files (chopping_and_ligation_data_*.npz)
#
# Outputs:
#   - Aggregated ligation maps
#   - Ligation probability vs genomic distance
#   - Ligation probability vs Euclidean distance
#   - Contact maps and distance matrices
#   - Ligation event statistics and time histories
#
# Paths and parameters should be adapted before execution.


#from joblib import Parallel, delayed
import numpy as np
import os.path
import scipy as sp
#import matplotlib.pyplot as plt
#import matplotlib as mpl
import warnings
warnings.filterwarnings('ignore')

nc = 5400
c = 7  # chromosome identity
r = '95.4-96.5'  # region identity
new_output_dir = '/path/to/aggregation_outputs'
output_dir   = '/path/to/ligation_outputs'
output_dir_2 = '/path/to/md_outputs'

nstart = 1
nfiles = 20
nsims = 1
gamma = 1

size = 5500     # number of beads in polymer model
rc = 1.5        # ligation cutoff distance
nsamples = 5000 # number of trajectories to include in ensemble
nframes = 5000  # max number of frames for time-history accumulation

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


ligation_prob_global_new = np.zeros((n_prob,dmax))
tot_nlig_new= np.zeros(n_prob)
all_fij_new = np.zeros((n_prob,len(x_edges)-1))
all_Dij_new = np.zeros((n_prob,len(x_edges)-1))
p_ligation_new = np.zeros((n_prob,len(x_edges)-1))             


ligation_map_global_new = np.zeros((n_prob,size, size))
ligation_distance_mat = []
for iprob in range(n_prob):
    ligation_distance_mat.append([])
net_eucledian_displacement_mat = []
for iprob in range(n_prob):
    net_eucledian_displacement_mat.append([])
ligation_time_mat = []
for iprob in range(n_prob):
    ligation_time_mat.append([])
cum_history_mat = []
for iprob in range(n_prob):
    cum_history_mat.append([])
accum_history = np.zeros((n_prob, nframes))

for filenumber in range(1,nfiles+1):
    for str_no in range(0,999,1):
    #for str_no in range(690,691,1):    
        for isim in range(1,nsims+1):
            fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(isim)
            ligationmap_file = output_dir+'/new_firsttime_ligation_map_'+fn_base+'_rc_'+str(rc)+'.npz'
            
            
            if (os.path.exists(ligationmap_file)) and (nf<nsamples):
                            
                print(f"Processing sample {nf+1}: {fn_base}")          

                data = np.load(ligationmap_file,allow_pickle=True)
                initial_distance_mat = data["initial_distance_mat"]
                initial_distance_final += initial_distance_mat
                init_contact_map_global = data["init_contact_map_global"]
                init_contact_map_final += init_contact_map_global
                
                #print(initial_distance_mat)
                
                
                #Calculation of original contact maps from original distance maps
                init_contact_map_global_2 = np.zeros((len(rthresh),size, size))
                for irt in range(n_rthresh):
                    distance_mat_thresholded        = rthresh[irt] - initial_distance_mat
                    init_contact_map_global_2[irt]  = (np.sign(distance_mat_thresholded)+1.)/2.
                init_contact_map_final_2 += init_contact_map_global_2

                
                g = np.load(output_dir_2+'/chopping_and_ligation_data_'+fn_base+'.npz')
                
                
                ligating_ends = g["ligating_ends"]
                ligating_ends = np.sort(np.append(ligating_ends,[0,size-1],0))
                n_sp = len(ligating_ends)
                n_spr = int(n_sp*(n_sp-1)/2)
            
            
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
                Log_book  = data["Log_book"]
                
                #Statistics are NOT capped by a number of ligation events. All ligations are considered. Each px will have different numbers of ligs.
                
                #Accumulate ligations and populate ligation maps up to the cap Ntotal_pmin
                for mm in range(n_prob):
               
            
                    Ntotal_px = int(Log_book[mm][-1][1])
                
                    #Populate ligation map with contacts up to the cap Ntotal_pmin
                    for lig_event in range(Ntotal_px):
                        
                        #Ligation ends
                        lig_end_1 = Log_book[mm][lig_event][2]
                        lig_end_2 = Log_book[mm][lig_event][3]
                        
                        #Add ligation event to ligation map
                        ligation_map_global_new[mm][lig_end_1][lig_end_2] += 1
        
                        #Add ligation event to statistics of initial Eucledian separation between ligating beads
                        initial_eucledian_separation = Log_book[mm][lig_event][4]
                        ligation_distance_mat[mm].append(initial_eucledian_separation)
                        #print(ligation_distance_mat[mm])
                        
                        #Add ligation event to statistics of total displacements of ligating beads until ligation
                        net_eucledian_displacement_1 = Log_book[mm][lig_event][5]
                        net_eucledian_displacement_2 = Log_book[mm][lig_event][6]
                        net_eucledian_displacement_mat[mm].append(net_eucledian_displacement_1)
                        net_eucledian_displacement_mat[mm].append(net_eucledian_displacement_2)
                        
                        #Add ligation event to statistics of ligation times
                        ligation_time = Log_book[mm][lig_event][0]
                        ligation_time_mat[mm].append(ligation_time)
                        
                        #Add ligation event to cumulative count of ligations vs time
                        accum_history[mm][ligation_time] +=1
                        
                        #print('mm', mm, 'Total', Ntotal_px, 'Lig No', lig_event + 1, 't', Log_book[mm][lig_event][0], 'ends', lig_end_1, lig_end_2, lig_end_1 < lig_end_2, ', initial euclid dist', Log_book[mm][lig_event][4], 'net disp.',  Log_book[mm][lig_event][5], Log_book[mm][lig_event][6])
                        

                    print('Total', np.sum(ligation_map_global_new[mm]), len(ligation_distance_mat[mm]), len(net_eucledian_displacement_mat[mm]), len(ligation_time_mat[mm]))

                    #print(ligation_distance_mat[mm])
                    #print(net_eucledian_displacement_mat[mm])
                    #print(ligation_time_mat[mm])
                    
                
                #Calculation of histograms WITHOUT CAP
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
    
#Calculation of history of ligation accumulations as a function of time
for mm in range(n_prob):
                    
    print('mm', mm)
    
    Ncumul = 0
    for tframe in range(nframes):
        ncount = accum_history[mm][tframe]
        if ncount != 0:
            Ncumul += ncount
            cum_history_mat[mm].append([tframe, Ncumul])    
    
    
#Calculation of histograms for each value of the probability with cap
for mm in range(n_prob):
                    
    print('mm', mm)
                    
    #Calculation of histogram of ligation probability vs. Eucledian distance.
    all_fij_new[mm] , x_edges = np.histogram(ligation_distance_mat[mm], bins = x_edges, density = False)
                    
    p_ligation_new[mm] = all_fij_new[mm]/all_Dij[mm]
    tot_nlig_new[mm] += np.sum(ligation_map_global_new[mm])    
    

init = 1
for mm in range(n_prob):
    ligation_map_global[mm] /= nf
    for i in range(dmax):
        ligation_prob_global[mm][i] =  np.mean(np.diagonal(ligation_map_global[mm], offset=(i+init)))
        
#Calculation of ligation prob. vs genomic distance        
init = 1
for mm in range(n_prob):
    ligation_map_global_new[mm] /= nf
    for i in range(dmax):
        ligation_prob_global_new[mm][i] =  np.mean(np.diagonal(ligation_map_global_new[mm], offset=(i+init)))

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


np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,init_contact_map_final_2=init_contact_map_final_2,rthresh=rthresh,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_L_final=Nlong_L_final,cycle_bins=cycle_bins,ligation_map_global_new=ligation_map_global_new,ligation_prob_global_new=ligation_prob_global_new,p_ligation_new=p_ligation_new,tot_nlig_new=tot_nlig_new,cum_history_mat=cum_history_mat)



#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,init_contact_map_final_2=init_contact_map_final_2,rthresh=rthresh,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_L_final=Nlong_L_final,cycle_bins=cycle_bins)

#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances, init_contact_map_final=init_contact_map_final,initial_distance_final=initial_distance_final,Pcycle_L_final=Pcycle_L_final,Ncycle_L_final=Ncycle_L_final,Nlong_final=Nlong_final,cycle_bins=cycle_bins)



#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=self_ligation_events,nf=nf,p_ligation1=p_ligation1,p_ligation2=p_ligation2,ligation_distances=ligation_distances)
#np.savez_compressed(mean_ligationmap_file,ligation_prob_global=ligation_prob_global,ligation_map_global=ligation_map_global,ligation_events=[tot_nlig,n_pairs],self_ligation_events=[tot_slig,n_lseg],length_ligation_2dhist=length_ligation_2dhist,dl=dl,dc=dc,nf=nf,h_length=h_length,ll_edges=ll_edges,ligation_distance_prob=ligation_distance_prob,ligation_distances=ligation_distances)
print(str(nf)+' files are used')
print(str(tot_nlig)+' ligation events out of '+str(n_pairs)+' open end pairs')
#print(str(tot_slig)+' total self ligation events out of '+str(n_lseg)+' long segments')


