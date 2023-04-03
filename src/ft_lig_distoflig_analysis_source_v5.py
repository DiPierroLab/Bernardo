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
lig_probs = np.array([1.0,0.1]) #,0.001,0.0001])
n_prob = len(lig_probs)

Ncounts_max = 1e15 #Maximum number of contacts allowed per experiment (for each ligation rate)

#Define threshold for max. number of successive time frames without ligations allowed after. If the threshold is exceeded, end ligating process
#for current value of p and move on to next probability.
Zero_ligs_thresh = 10000

input_dir = '@input_dir'
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
    
    
    #Start ligation experiment with smallest ligation rate and count total number of ligations,
    #which is not to exceed the maximum count Ncount.
    #Then run for higher ligation rates, stopping the experiment when the ligation count
    #matches the one obtained for the smallest rate.
    #In this way, the ligation maps will have the same total number of counts for all px.
    #Possible ligation events at a given time frame are visited in a random sequence to avoid bias.
    
    #Experiment for the smallest ligation rate.
    
    Ncounts_sm_pr    = 0 
    mm_smallest_prob = n_prob-1
    
    Nreject = 0
    ligation_map = np.zeros((n_prob,size,size))

    mm = mm_smallest_prob
    px = lig_probs[mm_smallest_prob]
    print('mm', mm, 'px', px)
    
    #Copy of available pairs list, to be decimated and updated as ligations take place
    ap_new = np.copy(ap[mm]) 

    #Keep track of the successive number of time frames during which there are no ligations registered.          
    Zero_ligs_reg = 0
    
    i = 0
    while((Ncounts_sm_pr < Ncounts_max) and (i<nt) and (Ncounts_sm_pr < n_seg) and (Zero_ligs_reg < Zero_ligs_thresh)) :
        s0 = gsd.hoomd.HOOMDTrajectory(f).read_frame(i)
        pos = s0.particles.position
        distance_mat = distance.cdist(pos, pos, 'euclidean')
        distance_mat_thresholded =  rc - distance_mat
        contact_map = (np.sign(distance_mat_thresholded)+1.)/2.
    
        contact_map_global += contact_map
        #ligation_map = np.zeros((n_prob,size,size))
        
        #Create a random order for the sequence of available pairs to remove bias
        ran_seq = list(range(int(n_avpr[mm])))
        random.shuffle(ran_seq)
        
        #Initialize number of ligations at time i
        N_lig_time_i = 0
            
        j = 0
        while (Ncounts_sm_pr < Ncounts_max) and (j < int(n_avpr[mm])):
            
            k = ran_seq[j] #Randomly selected available ligating-end pair
            
            if contact_map[tuple(ap[mm][k])] > 0.:
                
                z = np.random.random()
                if z <= px:
                    
                    #Check if pair of ends are already form a ligation
                    occupancy = ligation_map[mm][tuple(ap[mm][k])]
                    #print('occupancy', occupancy)
                    
                    #Check if any of the two ends have already ligated with any other end at all.
                    #If the first end is available because it hasn't ligated with any other end, then first_end_avail = 0. Likewise for the second end.
                    #If the first end has already ligated with some other end, then first_end_avail > 0. Likewise for the second end.
                    
                    first_end  = ap[mm][k][0]
                    second_end = ap[mm][k][1] 
                    
                    first_end_avail  = int(np.sum(ligation_map[mm][first_end][:])) + int(np.sum(ligation_map[mm][:][first_end]))
                    second_end_avail = int(np.sum(ligation_map[mm][second_end][:])) + int(np.sum(ligation_map[mm][:][second_end]))
                    
                    #In order for a ligation to happen, both ligating ends have to be completely available.                       
                    #No multi-ligations allowed for individual ligating ends.
                    
                    if (mm == mm_smallest_prob) and (Ncounts_sm_pr < Ncounts_max) and (occupancy == 0) and (first_end_avail == 0) and (second_end_avail == 0):
                            
                        #Cumulative number of ligations for current value of p    
                        Ncounts_sm_pr +=1
                        
                        #ligation_map[mm][tuple(ap[mm][k])] += 1
                        
                        ligation_map[mm][first_end][second_end] += 1
                        ligation_map[mm][second_end][first_end] += 1     
                        
                        #Count number of ligations in this time frame for this ligation probability
                        N_lig_time_i += 1 
                        
                        #Enter ligation into the log book
                        #Information entered: time frame i, ligation number, first bead, second bead, initial separation between beads,
                        #total displacement of first bead until ligation, total displacement of second bead until ligation 
                        
                        pos1i = initial_pos[first_end]
                        pos2i = initial_pos[second_end]
                        pos1f = pos[first_end]
                        pos2f = pos[second_end]
                        
                        displ1 = np.linalg.norm(pos1i-pos1f)
                        displ2 = np.linalg.norm(pos2i-pos2f)   
                        
                        Log_book[mm].append([i, Ncounts_sm_pr, first_end, second_end, initial_distance_mat[first_end, second_end], displ1, displ2]) 
                        
                        #print(Log_book[mm])
                        #input()
                        
                        
                        #Self-ligations / Cycles
                        
                        for m in range(n_lseg):
                            ancilla1 = [long_segments[m][0], long_segments[m][-1]]
                            ancilla2 = list(ap[mm][k])

                            ancilla1.sort()
                            ancilla2.sort()
                     
                            if ancilla1  == ancilla2:
                                self_ligation_events[mm][i][m] += 1
                                self_ligation_stats[mm][m] += len(long_segments[m])
       
                        #print('Ncounts_sm_pr', Ncounts_sm_pr)
        
        
                        #Scratching off available pairs whose ends match the ones of the ligation just effected
            
                        ap_new = np.delete(ap_new, np.where(ap_new==first_end)[0], axis=0)
                        ap_new = np.delete(ap_new, np.where(ap_new==second_end)[0], axis=0)                  
        
                    elif (occupancy > 0) or (first_end_avail > 0) or (second_end_avail > 0):
                        Nreject += 1
    
            j += 1
        
        #If there were ligations at time i, register total number of ligations that took place at i, and the cumulative number of ligs up to time i
        #Also, set the counter of the successive number of time steps without ligations to zero.
        if (N_lig_time_i > 0):
            N_lig_time[mm].append([i, N_lig_time_i, Ncounts_sm_pr])
            Zero_ligs_reg = 0
            
        #Keep track of number of time steps during which no ligations are recorded. If this number exceeds a threshold, the ligation process for probability p
        #ends and we proceed to the next probability.
        else:
            Zero_ligs_reg += 1
               
            #print(N_lig_time[mm])
            #input()
        
        #print('Zero ligs ', Zero_ligs_reg, 'i' , i)
        
        #Update the list of available ligating pairs by removing ligated pairs
        #This is done to reflect the fact that ligated ends are not available for any further ligations in the future.
        #print('preupdate',  ap[mm]==ap_new, len(ap[mm]))
        ap[mm] = np.copy(ap_new)

        n_avpr[mm] = len(ap[mm])
        
        #print('postupdate', ap[mm]==ap_new, len(ap[mm]))        
    
        Ntotal += 1
        if i % 100 == 0:
            print("Reading frame {:} of {:}".format(i, f.nframes), 'Cum Number of ligations', Ncounts_sm_pr, 'Cum Number of rejections', Nreject, 'Counts in lig map', np.sum(ligation_map[mm]), 'Max entry in map', np.max(ligation_map[mm]), 'Max no. neighbors', max(np.sum(ligation_map[mm],axis=0)), 'Number of segments', n_seg)
            #print(Log_book[mm])
            #input()
            
        i += 1 
        
    #Update ligation map  for this value of ligation rate. We keep only the upper triangle of the symmetric ligation map.
    ligation_map_with_prob[mm] += np.triu(ligation_map[mm]) 
    
    #Record total number of time frames elapsed for this value of ligation rate.
    Max_lig_time[mm] = i
        
    print('mm ', mm ,'px ', lig_probs[mm], 'Sum of ligation map ', np.sum(ligation_map_with_prob[mm]), 'Number of ligs', Ncounts_sm_pr, 'Number of rejects', Nreject, 'Counts in final lig map', np.sum(ligation_map_with_prob[mm]), 'Max entry in final map', np.max(ligation_map_with_prob[mm]), 'Max no. neighbors in final map', max(np.sum(ligation_map_with_prob[mm],axis=0)), 'Number of segments', n_seg)
    
    #Experiments for the remaining (higher) probability rates of ligation
    #Maximum number of ligations to count is capped by the total number of ligations measured for the smallest px.
    
    Ncounts_sm_pr = Ncounts_max #10000000
    
    if(n_prob > 1):
    
        for mm in range(n_prob-2,-1,-1):
            px = lig_probs[mm]
            Ncount = 0
            Nreject = 0
            print('mm', mm, 'px', px)
            
            #Copy of available pairs list, to be decimated and updated as ligations take place
            ap_new = np.copy(ap[mm]) 
            
            #Keep track of the successive number of time frames during which there are no ligations registered.          
            Zero_ligs_reg = 0
    
            i = 0
 
            while((Ncount < Ncounts_sm_pr) and (i<nt) and (Ncount < n_seg) and (Zero_ligs_reg < Zero_ligs_thresh)) :
                s0 = gsd.hoomd.HOOMDTrajectory(f).read_frame(i)
                pos = s0.particles.position
                distance_mat = distance.cdist(pos, pos, 'euclidean')
                distance_mat_thresholded =  rc - distance_mat
                contact_map = (np.sign(distance_mat_thresholded)+1.)/2.
    
                contact_map_global += contact_map
                #ligation_map = np.zeros((n_prob,size,size))
            
                #Create a random order for the sequence of available pairs to remove bias
                ran_seq = list(range(int(n_avpr[mm])))
                random.shuffle(ran_seq) 
                
                #Initialize number of ligations at time i
                N_lig_time_i = 0
            
                j = 0
                while (Ncount < Ncounts_sm_pr) and (j < int(n_avpr[mm])):
                
                    k = ran_seq[j] #Randomly selected available ligating-end pair
                
                    if contact_map[tuple(ap[mm][k])] > 0.:

                        z = np.random.random()
                        if z <= px:
                        
                            #Check if pair of ends are already form a ligation
                            occupancy = ligation_map[mm][tuple(ap[mm][k])]
                    
                            #Check if any of the two ends have already ligated with any other end at all.
                            #If the first end is available because it hasn't ligated with any other end, then first_end_avail = 0. Likewise for the second end.
                            #If the first end has already ligated with some other end, then first_end_avail > 0. Likewise for the second end.
                    
                            first_end  = ap[mm][k][0]
                            second_end = ap[mm][k][1] 
                    
                            first_end_avail  = np.sum(ligation_map[mm][first_end][:]) + np.sum(ligation_map[mm][:][first_end])
                            second_end_avail = np.sum(ligation_map[mm][second_end][:]) + np.sum(ligation_map[mm][:][second_end]) 
                    
                            #In order for a ligation to happen, both ligating ends have to be completely available.                       
                            #No multi-ligations allowed for individual ligating ends.
    
                            if (mm != mm_smallest_prob) and (Ncount < Ncounts_sm_pr) and (occupancy == 0) and (first_end_avail == 0) and (second_end_avail == 0):
                                
                                Ncount +=1
                            
                                ligation_map[mm][first_end][second_end] += 1
                                ligation_map[mm][second_end][first_end] += 1     
                                
                                #Count number of ligations in this time frame for this ligation probability
                                N_lig_time_i += 1 
                                
                                #Enter ligation into the log book
                                #Information entered: time frame i, ligation number, first bead, second bead, initial separation between beads,
                                #total displacement of first bead until ligation, total displacement of second bead until ligation 
                        
                                pos1i = initial_pos[first_end]
                                pos2i = initial_pos[second_end]
                                pos1f = pos[first_end]
                                pos2f = pos[second_end]
                        
                                displ1 = np.linalg.norm(pos1i-pos1f)
                                displ2 = np.linalg.norm(pos2i-pos2f)   
                        
                                Log_book[mm].append([i, Ncounts_sm_pr, first_end, second_end, initial_distance_mat[first_end, second_end], displ1, displ2]) 
                        
                                #print(Log_book[mm])
                                #input()
                                
                                #Self-ligations / Cycles
                                
                                for m in range(n_lseg):
                                    ancilla1 = [long_segments[m][0], long_segments[m][-1]]
                                    ancilla2 = list(ap[mm][k])

                                    ancilla1.sort()
                                    ancilla2.sort()
                     
                                    if ancilla1  == ancilla2:
                                        self_ligation_events[mm][i][m] += 1
                                        self_ligation_stats[mm][m] += len(long_segments[m])

                                #print('Ncount', Ncount, 'Ncounts_sm_pr', Ncounts_sm_pr)
                                
                                
                                #Scratching off available pairs whose ends match the ones of the ligation just effected
            
                                ap_new = np.delete(ap_new, np.where(ap_new==first_end)[0], axis=0)
                                ap_new = np.delete(ap_new, np.where(ap_new==second_end)[0], axis=0)                            
                            
                            elif (occupancy > 0) or (first_end_avail > 0) or (second_end_avail > 0):
                                Nreject += 1 
                            
                    j +=1

                    
                #If there were ligations at time i, register total number of ligations that took place at i, and the cumulative number of ligs up to time i
                #Also, set the counter of the successive number of time steps without ligations to zero.
                if (N_lig_time_i > 0):
                    N_lig_time[mm].append([i, N_lig_time_i, Ncount])    
                    Zero_ligs_reg = 0
            
                #Keep track of number of time steps during which no ligations are recorded. If this number exceeds a threshold, the ligation process for probability p ends and we proceed to the next probability.
                else:
                    Zero_ligs_reg += 1
               
                    #print(N_lig_time[mm])
                    #input()
        
                #print('Zero ligs ', Zero_ligs_reg, 'i' , i)
                    
                #Update the list of available ligating pairs by removing ligated pairs
                #This is done to reflect the fact that ligated ends are not available for any further ligations in the future.
                #print('preupdate',  ap[mm]==ap_new, len(ap[mm]))
                ap[mm] = np.copy(ap_new)

                n_avpr[mm] = len(ap[mm])
                #print('postupdate', ap[mm]==ap_new, len(ap[mm]))        
                    
                Ntotal += 1
                if i % 100 == 0:
                    #print("Reading frame {:} of {:}".format(i, f.nframes))         
                    print("Reading frame {:} of {:}".format(i, f.nframes), 'Cum Number of ligations', Ncount, 'Cum Number of rejections', Nreject, 'Counts in lig map', np.sum(ligation_map[mm]), 'Max entry in map', np.max(ligation_map[mm]), 'Max no. neighbors', max(np.sum(ligation_map[mm],axis=0)), 'Number of segments', n_seg)      
            
                i += 1
    
            #Update ligation map  for this value of ligation rate 
            ligation_map_with_prob[mm] += np.triu(ligation_map[mm])
            
            #Record total number of time frames elapsed for this value of ligation rate.
            Max_lig_time[mm] = i
        
            print('mm ', mm ,'px ', lig_probs[mm], 'Sum of ligation map ', np.sum(ligation_map_with_prob[mm]), 'Number of ligs', Ncounts_sm_pr, 'Number of rejects', Nreject, 'Counts in final lig map', np.sum(ligation_map_with_prob[mm]), 'Max entry in final map', np.max(ligation_map_with_prob[mm]), 'Max no. neighbors in final map', max(np.sum(ligation_map_with_prob[mm],axis=0)),  'Number of segments', n_seg)        
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
