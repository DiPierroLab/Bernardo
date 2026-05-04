#!/usr/bin/env python
# coding: utf-8

# Example generated script produced from ft_lig_distoflig_analysis_source_v6.py
# using gen_firsttime_ligation_maps_from_manystr_v6.sh.
#
# This file illustrates the structure of a generated ligation-analysis script
# for one MD trajectory. It is included for documentation only; paths should
# be adapted before execution.

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
import copy
warnings.filterwarnings('ignore')

seed_num = 29554
random.seed(seed_num)
c = 7  # chromosome identity
r = '95.4-96.5'  # region identity
filenumber = 1
str_no = 0
run_no = 1
output_dir = '/path/to/md_outputs'
new_output_dir = '/path/to/ligation_analysis_outputs'
fn_base = 'C_7_region_95.4-96.5_nkeep_5400_str_0_n_1_sim_1'
rc = 1.5
lig_probs = np.array([1.0,0.1,0.01,0.001])
n_prob = len(lig_probs)

Ncounts_max = 1e15 #Maximum number of contacts allowed per experiment (for each ligation rate)

#Define threshold for max. number of successive time frames without ligations allowed after. If the threshold is exceeded, end ligating process
#for current value of p and move on to next probability.
Zero_ligs_thresh = 3000

input_dir = '/path/to/minimized_structures'
minimized_str_file = input_dir+'/minimized_C_'+str(c)+'_region_'+r+'_file_'+str(filenumber)+'_str_'+str(str_no)+'_minimized.gsd'
f0 = gsd.fl.open(minimized_str_file,'rb')
initial_config = gsd.hoomd.HOOMDTrajectory(f0).read_frame(f0.nframes-1)
pos = initial_config.particles.position
initial_distance_mat = distance.cdist(pos, pos, 'euclidean')

trajfile = output_dir+'/chromosome_crosslinked_fixed_chopped_'+fn_base+'.gsd'
contactmap_file = new_output_dir+'/new_firsttime_ligation_map_'+fn_base+'_rc_'+str(rc)+'.npz'

#Initial contact map, over non-chopped, non-ligated original structure

size = initial_config.particles.N
init_contact_map_global  = np.zeros((size, size))
distance_mat_thresholded = rc - initial_distance_mat
init_contact_map_global  = (np.sign(distance_mat_thresholded)+1.)/2.

#Initial positions of the particles for future calculation of displacements
initial_pos = initial_config.particles.position

#Ligation process

if not os.path.exists(contactmap_file):
    f = gsd.fl.open(trajfile,'rb')
    s0 = gsd.hoomd.HOOMDTrajectory(f).read_frame(0)
    nt = f.nframes
    size = s0.particles.N
        
    g = np.load(output_dir+'/chopping_and_ligation_data_'+fn_base+'.npz')
    ligating_ends = g["ligating_ends"]
    avail_ends    = np.copy(ligating_ends)
    avail_ends    = np.append(avail_ends, size-1)
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
    
    #Add final segment        
    e1 = rbs[len(rbs)-1]
    seg1 = np.arange(e1, size)
    segments.append(seg1)
    
    n_seg = len(segments)
    n_lseg = len(long_segments)
    
    #Maximum number of contacts allowed per experiment (for each ligation rate).
    #Set it equal to number of segments, because there can't be more ligations than this number if there are no multiligations.
    Ncounts_max =  n_seg  #1e15 
    
    contact_map_global = np.zeros((size, size))
    #ligation_map_global = np.zeros((size, size))
    ligation_map_with_prob = np.zeros((n_prob,size,size))
    self_ligation_events = np.zeros((n_prob,nt,n_lseg))
    
    self_ligation_stats  = np.zeros((n_prob,n_lseg))
        
    #List that contains the 2 ends of each long segment, useful for assessing cycles
    long_segment_ends = []
    for i in range(n_lseg):
        long_segment_ends.append([long_segments[i][0], long_segments[i][-1]])
        
    Ntotal = 0
    ligated_ends = []
    available_pairs = ligating_pairs
    ap = []
    n_avpr = np.zeros(n_prob)
    for i in range(n_prob):
        ap.append(np.array(available_pairs))
        n_avpr[i] = len(available_pairs)        
    
    avends = []
    n_avends = np.zeros(n_prob)
    for i in range(n_prob):
        avends.append(np.array(avail_ends))
        n_avends[i] = len(avail_ends)

    dx = 0.2
    x_edges = np.arange(dx,10+dx,dx)
    fij_global = np.zeros((n_prob,len(x_edges)-1))
    x_centers = 0.5*(x_edges[:-1]+x_edges[1:])    

    #Create a list to store the history of the number of ligations as a function of time for each ligation prob.
    #This is a list with n_prob lists, one for each lig. probability.
    #For each prob, the time frame, the number of ligations in that time frame and the cumulative number of ligations up to that time frame are recorded
    N_lig_time = []
    for iprob in range(n_prob):
        N_lig_time.append([])
        
    #Create a list to store a log records of each ligation that took place, including time at which ligation happened, the beads that ligated,
    #the original separation (Euclidean distance) and the total distance traveled by each bead until ligation.
    #This is a list with n_prob lists, one for each lig. probability.
    Log_book= []
    for iprob in range(n_prob):
        Log_book.append([])
    
    #Create a list to store the final time frame reached for each ligation probability. This is because the ligation process for a given lig. prob.
    #will be killed after a certain threshold of successive frames without ligations (dormant) is exceeded.
    Max_lig_time = np.zeros(n_prob)
    
    #Number of rejections for different ligation rates
    Nreject = np.zeros(n_prob)
    
    #List of ligations to replace adjacency matrix. This list contains as many entries as there are beads.
    #For each bead, it contains the label of the bead to which it has ligated. If the bead has not ligated yet, then its entry is False.
    lig_list = []
    for iprob in range(n_prob):
        lig_list.append([False for _ in range(size)])
        
    #Keep track of the successive number of time frames during which there are no ligations registered.          
    Zero_ligs_reg = np.zeros(n_prob)
    
    #Initialize ligation maps for different ligation rates
    ligation_map = np.zeros((n_prob,size,size))
    
    #Initialize cumulative number of counts
    Ncount = np.zeros(n_prob)
    
    #Number of available pairs at each time frame
    n_avail_pairs = np.zeros(n_prob)            

    #Run ligation experiment for each time frame, exploring each of the ligating rates.
    #Stop process if the number of ligations is zero after a number consecutive frames equal to Zero_ligs_thresh
    
    i = 10

    while( (i<nt) and (min(Ncount) < Ncounts_max) and (min(Ncount) < n_seg) and (min(Zero_ligs_reg) < Zero_ligs_thresh)) :
        
        s0 = gsd.hoomd.HOOMDTrajectory(f).read_frame(i)
        pos = s0.particles.position
        distance_mat = distance.cdist(pos, pos, 'euclidean')
        distance_mat_thresholded =  rc - distance_mat
        contact_map = (np.sign(distance_mat_thresholded)+1.)/2.
        
        #contact_map_global += contact_map
        
        #Copy of available pairs list for each prob, to be decimated and updated as ligations take place
        ap_new = copy.deepcopy(ap)
        
        #print(ap_new)
        #print(ap)
    
        #Copy of available ligating ends list, to be decimated and updated as ligations take place
        avail_ends_new   = copy.deepcopy(avends)
        n_avail_ends_new = copy.deepcopy(n_avends)
        
        #print(avail_ends_new)
        #print(avends)
        
        #Initialize number of ligations at time i
        N_lig_time_i = np.zeros(n_prob)
        
        for mm in range(n_prob):
            
            if((Ncount[mm] < Ncounts_max) and (Ncount[mm] < n_seg) and (Zero_ligs_reg[mm] < Zero_ligs_thresh)) :
        
                px = lig_probs[mm]
                #print('i', i, 'mm', mm, 'px', px)
    
                #Create a list that contains only the pairs of beads that are within physical proximity at this point in time.
                #Do this extracting only the non-zero entries of the upper triangular matrix (without the diagonal) of the contact map to avoid overcounting.
                #Set to zero all columns and rows that correspond to beads that are NOT in the list of available ends.
                #From this matrix, a list is compiled that only counts beads that are available ends within proximity at this point in time. 
        
                tria_cont = np.triu(contact_map-np.diag(np.diag(contact_map)))  #Count only pairs of beads in proximity    
                for ibead in range(size): #Filter out beads that are not segment ends
                    if ((ibead in avail_ends_new[mm]) == False):
                        tria_cont[ibead]   = np.zeros(size)
                        tria_cont[:,ibead] = np.zeros(size)          
                phys_prox_list    = np.transpose(np.nonzero(tria_cont))
                n_avail_pairs[mm] = len(phys_prox_list) 

                #print(n_avail_pairs[mm])                
                #print(phys_prox_list)

                #Create a random order for the list of available pairs at this point in time to remove bias.
                rand_seq = list(range(int(n_avail_pairs[mm])))        
                random.shuffle(rand_seq)

                #Loop over available pairs at time frame i
                for j in range(int(n_avail_pairs[mm])):

                    k = rand_seq[j] #Randomly selected available ligating-end pair
                    
                    #print(j) #print(k) #print(phys_prox_list[k]) #print(tuple(phys_prox_list[k])) #print(contact_map[tuple(phys_prox_list[k])]) #sys.exit()
            
                    if contact_map[tuple(phys_prox_list[k])] > 0.:                     

                    #print('post', j, k)
                
                        z = np.random.random()
                        if z <= px:
                    
                            #Check if any of the two ends have already ligated with any other end at all.
                            #If the first end is available because it hasn't ligated with any other end, then first_end_avail = 0. Likewise for the second end.
                            #If the first end has already ligated with some other end, then first_end_avail > 0. Likewise for the second end.
                    
                            #first_end  = ap[mm][k][0] #second_end = ap[mm][k][1] 
                    
                            first_end  = phys_prox_list[k][0]
                            second_end = phys_prox_list[k][1] 
                                        
                            #print('end labels', first_end, second_end)
                    
                            #first_end_avail  = int(np.sum(ligation_map[mm][first_end][:])) + int(np.sum(ligation_map[mm][:][first_end]))
                            #second_end_avail = int(np.sum(ligation_map[mm][second_end][:])) + int(np.sum(ligation_map[mm][:][second_end]))
                    
                            #Check if pair of ends are already form a ligation
                            #occupancy = ligation_map[mm][tuple(ap[mm][k])]
                            occupancy = ligation_map[mm][tuple(phys_prox_list[k])]
                            occ_bead1 = lig_list[mm][first_end] 
                            occ_bead2 = lig_list[mm][second_end]
                    
                            #print('liglistlength', len(lig_list)) #print('occupancy', occupancy, occ_bead1, occ_bead2, first_end, second_end )
                            #print(mm==mm_smallest_prob, Ncounts_sm_pr < Ncounts_max, occupancy ==0, first_end_avail == 0, second_end_avail == 0)
                    
                            #In order for a ligation to happen, both ligating ends have to be completely available.                       
                            #No multi-ligations allowed for individual ligating ends.  

                            if (Ncount[mm] < Ncounts_max) and (occupancy == 0) and (occ_bead1 == False) and (occ_bead2 == False):
                        
                                #Cumulative number of ligations for current value of p    
                                Ncount[mm] +=1
                        
                                #ligation_map[mm][tuple(ap[mm][k])] += 1 #print('lig map pre', ligation_map[mm][first_end][second_end], ligation_map[mm][second_end][first_end] )
                        
                                ligation_map[mm][first_end][second_end] += 1
                                ligation_map[mm][second_end][first_end] += 1   
                        
                                #print('lig map post', ligation_map[mm][first_end][second_end], ligation_map[mm][second_end][first_end] ) #print('preupdate', lig_list[mm][first_end], lig_list[mm][second_end], first_end, second_end)
                        
                                lig_list[mm][first_end]  = second_end
                                lig_list[mm][second_end] = first_end
                        
                                #print('postupdate', lig_list[mm][first_end], lig_list[mm][second_end], first_end, second_end)
                            
                        
                                #Count number of ligations in this time frame for this ligation probability
                                N_lig_time_i[mm] += 1 
                                #print('N_lig_time_i', N_lig_time_i[mm])
                        
                                #Enter ligation into the log book
                                #Information entered: time frame i, ligation number, first bead, second bead, initial separation between beads,
                                #total displacement of first bead until ligation, total displacement of second bead until ligation 
                        
                                pos1i = initial_pos[first_end]
                                pos2i = initial_pos[second_end]
                                pos1f = pos[first_end]
                                pos2f = pos[second_end]
                        
                                displ1 = np.linalg.norm(pos1i-pos1f)
                                displ2 = np.linalg.norm(pos2i-pos2f)   
                        
                                Log_book[mm].append([i, Ncount[mm], first_end, second_end, initial_distance_mat[first_end, second_end], displ1, displ2]) 
                        
                                #print(Log_book[mm][-1]) #sys.exit() #input()
                        
                                #Self-ligations / Cycles
                                ancilla = [first_end, second_end]
                                ancilla.sort()
                                
                                if ancilla in long_segment_ends:
                                    ancilla_index = long_segment_ends.index(ancilla)
                                    self_ligation_events[mm][i][ancilla_index] += 1
                                    self_ligation_stats[mm][ancilla_index] += len(long_segments[ancilla_index])
        
                                #Scratching off available pairs whose ends match the ones of the ligation just effected
            
                                ap_new[mm] = np.delete(ap_new[mm], np.where(ap_new[mm]==first_end)[0], axis=0)
                                ap_new[mm] = np.delete(ap_new[mm], np.where(ap_new[mm]==second_end)[0], axis=0)   
                                #print(ap_new[mm]) #print(len(ap_new[mm]))
                    
                                avail_ends_new[mm] = np.delete(avail_ends_new[mm], np.where(avail_ends_new[mm]==first_end)[0], axis=0)
                                avail_ends_new[mm] = np.delete(avail_ends_new[mm], np.where(avail_ends_new[mm]==second_end)[0], axis=0)
                                #print(avail_ends_new[mm]) #print(len(avail_ends_new[mm]))
        
                            elif (occupancy > 0) or (occ_bead1 != False) or (occ_bead2 != False):
                                Nreject[mm] += 1
                                #print('Nreject[mm]',Nreject[mm]) #print('________________________________________')

                #If there were ligations at time i, register total number of ligations that took place at i, and the cumulative number of ligs up to time i
                #Also, set the counter of the successive number of time steps without ligations to zero.
                if (N_lig_time_i[mm] > 0):
                    N_lig_time[mm].append([i, N_lig_time_i[mm], Ncount[mm]])
                    Zero_ligs_reg[mm] = 0
                    
                    #print('N_lig_time_i[mm]', N_lig_time_i[mm]) #print('Zero_ligs_reg[mm]', Zero_ligs_reg[mm])
            
                #Keep track of number of time steps during which no ligations are recorded. If this number exceeds a threshold, the ligation process for probability p
                #ends and we proceed to the next probability.
                else:
                    Zero_ligs_reg[mm] += 1
                    
                    #print('N_lig_time_i[mm]', N_lig_time_i[mm]) #print('Zero_ligs_reg[mm]', Zero_ligs_reg[mm]) #print(N_lig_time[mm]) #input()
                #print('Zero ligs ', Zero_ligs_reg, 'i' , i)
        
                #Update the list of available ligating pairs by removing ligated pairs
                #This is done to reflect the fact that ligated ends are not available for any further ligations in the future.
                #print('preupdate',  ap[mm]==ap_new, len(ap[mm]))
                ap[mm]       = np.copy(ap_new[mm])
                n_avpr[mm]   = len(ap[mm])
        
                avends[mm]   = np.copy(avail_ends_new[mm])
                n_avends[mm] = len(avends[mm])
        
                #print('available ends post', n_avends[mm]) #print('available ends pre',  n_avail_ends_new) #print('postupdate', ap[mm]==ap_new, len(ap[mm]))      
    
                #Record total number of time frames elapsed for this value of ligation rate.
                Max_lig_time[mm] = i
                #print(avail_ends_new)
                        
        if i % 200 == 0:
            for mm in range(n_prob):
                print("Reading frame {:} of {:}".format(i, f.nframes), 'mm', mm, 'px', lig_probs[mm], 'Cum Number of ligations', Ncount[mm], 'Cum Number of rejections', Nreject[mm], 'Counts in lig map', np.sum(ligation_map[mm]), 'Max entry in map', np.max(ligation_map[mm]), 'Max no. neighbors', max(np.sum(ligation_map[mm],axis=0)), 'Number of segments', n_seg, 'Init No. of avail. ends: ', n_avail_ends_new[mm], 'Current No. of avail. ends: ', n_avends[mm])
            #print(Log_book[mm])
            #input()
    
        Ntotal += 1
        i += 1 
        
    #Update ligation map  for this value of ligation rate. We keep only the upper triangle of the symmetric ligation map.
    ligation_map_with_prob += np.triu(ligation_map) 

    for mm in range(n_prob):
        print('i', i, 'mm ', mm ,'px ', lig_probs[mm], 'Sum of ligation map ', np.sum(ligation_map_with_prob[mm]), 'Number of ligs', Ncount[mm], 'Number of rejects', Nreject[mm], 'Counts in final lig map', np.sum(ligation_map_with_prob[mm]), 'Max entry in final map', np.max(ligation_map_with_prob[mm]), 'Max no. neighbors in final map', max(np.sum(ligation_map_with_prob[mm],axis=0)), 'Number of segments', n_seg, 'Init No. of avail. ends: ', n_avail_ends_new[mm], 'Final No. of avail. ends: ', n_avends[mm])
    
    for mm in range(n_prob):
        print('mm ', mm ,'px ', lig_probs[mm], 'Sum of ligation map ', np.sum(ligation_map_with_prob[mm]), 'Max. lig. time ', Max_lig_time[mm])                        
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
        
        
    #Contact probabilities of original structure vs. genomic distance and Eucledian distance
    contact_distance_mat = np.multiply(init_contact_map_global, distance_mat_thresholded)
    fij_o, x_edges       = np.histogram(contact_distance_mat.flat, bins=x_edges, density=False)
    
    #Self-ligation / cyclization statistics
    #& Statistics of lengths of long segments
    
    cycle_bins = np.arange(5,101) - 0.5
    
    len_long_seg = np.zeros(n_lseg)
    for i in range(n_lseg):
        len_long_seg[i] = len(long_segments[i])
        
    Nlong_L, cycle_bins = np.histogram(len_long_seg, bins=cycle_bins, density=False)
    
    Pcycle_L = np.zeros((n_prob, len(cycle_bins)-1))
    for i in range(n_prob):
        cycle_bins = np.arange(5,101) - 0.5
        Pcycle_L[i], cycle_bins = np.histogram(self_ligation_stats[i], bins=cycle_bins, density=False)
           
np.savez_compressed(contactmap_file,initial_distance_mat=initial_distance_mat,init_contact_map_global=init_contact_map_global,ligation_map=ligation_map_with_prob,Ntotal=Ntotal,self_ligation_events=self_lig_ev,long_segments=long_segments,Dij=Dij,fij_global=fij_global,fij_o=fij_o,contact_distance_mat=contact_distance_mat,x_edges=x_edges,Nlong_L=Nlong_L,Pcycle_L=Pcycle_L,cycle_bins=cycle_bins,N_lig_time=N_lig_time,Log_book=Log_book,Max_lig_time=Max_lig_time)                               
     
#np.savez_compressed(contactmap_file,initial_distance_mat=initial_distance_mat,init_contact_map_global=init_contact_map_global,ligation_map=ligation_map_with_prob,Ntotal=Ntotal,self_ligation_events=self_lig_ev,long_segments=long_segments,Dij=Dij,fij_global=fij_global,x_edges=x_edges,Nlong_L=Nlong_L,Pcycle_L=Pcycle_L,cycle_bins=cycle_bins) 

#np.savez_compressed(contactmap_file,ligation_map=ligation_map_with_prob,Ntotal=Ntotal,self_ligation_events=self_lig_ev,long_segments=long_segments,Dij=Dij,fij_global=fij_global,x_edges=x_edges)
#    
    
    
    
    
    
    
    
    
    
    
 
