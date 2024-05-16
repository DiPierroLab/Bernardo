#!/usr/bin/env python
# coding: utf-8

# In[1]:


import hoomd, hoomd.md
import mdtraj as md
import mbuild as mb
import numpy  as np
import random
#from joblib import Parallel, delayed
import multiprocessing
import os.path
import warnings
warnings.filterwarnings('ignore')
import sys
import gsd
import gsd.fl, gsd.hoomd
import math
import copy

#new FENE potential class
def fene2(r, rmin, rmax, kb, r0, rshift):
    #print(r)

    if r < rshift or r >= r0 + rshift: #a == 1.0:
        V = 0.0
        F = 0.0
    else:
        a = (r - rshift)/r0
        x = 1.0 - a*a
        V = -0.5*kb*r0*r0*np.log(x)
        F = -kb*r0*a/x

    return (V, F)

#def fenewca(r, rmin, rmax, kb, r0, rshift, epsilon, sigma):
##     print(r)
#    a = (r - rshift)/r0
#    x = 1.0 - a*a
#    b = sigma/r
#    r6 = b**6.0
#    r12 = r6*r6
#    r7 = b**7.0
#    r13 = r6*r7
#
#    V = -0.5*kb*r0*r0*np.log(x) + 4.0*epsilon*(r12 - r6 + 1.0/4.0)
#    if a == 1.0:
#        F = 0.0
#    else:
#        F = -kb*r0*a/x + (48.0*epsilon/sigma)*(r13 - 0.5*r7)
#
#    return (V, F)

#reading input file containing output of Langevin dynamics simulation with energy minimizaed polymer 
hoomd.context.initialize("");
minimized_str_file = '@minimized_str_file'
f = gsd.fl.open(minimized_str_file,'rb')

#parameters
ncut = @ncut  #number of bonds to cut after each time evolution - control degree of degradation
nsteps = @nsteps #number of time steps to evolve the structure between cuts
ncycles = @ncycles #number of cycles of evolution and degradation
c = @chr  # chromosome identity
r = '@region'  # region identity
filenumber = @filenumber
str_no = @str
run_no = @run
output_dir = '@output_dir'
#fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(run_no)

#creating hoomd snapshot with final frame
snapshot = hoomd.data.gsd_snapshot(minimized_str_file,frame=f.nframes-1)
n_particles = snapshot.particles.N

random.seed(@seed1)

#nc = nkeep # number of bonds to keep after chopping
#all_bond_list = np.arange(n_particles-1)
#bond_list = np.sort(np.random.choice(range(n_particles-1), nc, replace=False))
#removed_bond_list = np.setdiff1d(all_bond_list,bond_list)
#nr = len(removed_bond_list)
#snapshot.bonds.resize(nc)

#Initial number of bonds in the polymer
nc = (n_particles-1)

#Initial list with all bonds
init_bond_list = range(n_particles-1)
init_all_bond_list = np.arange(n_particles-1)

#base for file names
fn_base_static = 'C_'+str(c)+'_region_'+str(r)+'_ncut_'+str(ncut)+'_nsteps_'+str(nsteps)+'_ncycle_'+str(ncycles)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(run_no)

#Arrays that will store the special particles where things are chopped, the bonds kept and the angles kept at each cycle
special_particles_record = []
kept_bonds_record        = []
kept_angles_record       = []

#Arrays that will store the removed angle and bond lists
removed_angles_record = []
removed_bonds_record  = []
fin_bonds_record      = []

for icycle in range(ncycles):
    
    #number of bonds kept after chopping at this stage
    nc = nc - ncut
    
    fin_bond_list = np.sort(np.random.choice(init_bond_list, nc,replace=False)) #list of bonds kept after chopping
    removed_bond_list = np.setdiff1d(init_all_bond_list,fin_bond_list) #list of bonds removed after chopping
    nr = len(removed_bond_list) #number of removed bonds
    snapshot.bonds.resize(nc)

    #updating bonds list
    all_bonds = []
    i1 = 0
    for i in range(n_particles-1):
        i2 = i1+1
        all_bonds.append([i1,i2])
        i1 = i2

    for i in range(nc):
        snapshot.bonds.group[i] = all_bonds[fin_bond_list[i]]
    
    special_particles = []
    for i in range(len(removed_bond_list)):
        m = removed_bond_list[i]
        special_particles.append(all_bonds[m][0])
        special_particles.append(all_bonds[m][1])    

    #getting the identities of all beads whose bonds are being chopped - these are the open-ends
    special_particles = np.unique(special_particles)
    n_sp = len(special_particles)
    
    #preparing the list of removed angles and updating the angles list in hoomd snapshot
    removed_angles = []
    for i in range(nr):
        bond_id = removed_bond_list[i]
        if bond_id ==0:
            removed_angles.append(0)
        elif bond_id ==n_particles-1:
            removed_angles.append(n_particles-2)
        else:
            removed_angles.append(bond_id-1)
            removed_angles.append(bond_id)

    all_angle_list = np.arange(n_particles-2)
    angles_list = np.setdiff1d(all_angle_list,removed_angles)
    n_angles = len(angles_list)

    all_angles = []
    i1 = 0
    for i in range(n_particles-2):
        i2 = i1+1
        i3 = i2+1
        all_angles.append([i1,i2,i3])
        i1 = i2
        
    snapshot.angles.resize(n_angles)

    for i in range(n_angles):
        snapshot.angles.group[i] = all_angles[angles_list[i]]

    init_bond_list = copy.deepcopy(fin_bond_list)
    
    #Storing list of free ends, kept bonds and kept angles in this cycle
    special_particles_record.append(copy.deepcopy(special_particles))
    kept_bonds_record.append(copy.deepcopy(snapshot.bonds.group))
    kept_angles_record.append(copy.deepcopy(snapshot.angles.group))
    removed_bonds_record.append(copy.deepcopy(removed_bond_list))
    removed_angles_record.append(copy.deepcopy(removed_angles))
    fin_bonds_record.append(copy.deepcopy(fin_bond_list))
    
    #using the updated snapshot to initialize the system to be simulated    
    if icycle == 0:
        system = hoomd.init.read_snapshot(snapshot)
    else:
        system.restore_snapshot(snapshot)

    #declaring the neighbor list using cell list method
    if icycle == 0:
        nl = hoomd.md.nlist.cell();
        nl.reset_exclusions(exclusions = []);

        #declaring the new fene potential using the table method of hoomd.md.pair
        fene = hoomd.md.bond.table(width=1000)
        fene.bond_coeff.set('polymer',func=fene2, rmin=0., rmax=5.0, coeff = dict(kb=30, r0=1.5, rshift=1.5))

        #declaring the harmonic angle potential
        #angle = hoomd.md.angle.harmonic()
        #angle.angle_coeff.set('backbone',k=4,t0=np.pi)

        #declaring the trucanted Lennard-Jones potential (Weeks, Chandler and Anderson (WCA) version)
        lj = hoomd.md.pair.lj(r_cut=3.0, nlist=nl)
        lj.pair_coeff.set('CA', 'CA', epsilon=1.0, sigma=1.0, r_cut=2**(1.0/6.0));

        #group containing all beads
        all = hoomd.group.all()

        #group containing all evolving beads
        evolving_tag_list = np.arange(n_particles)
        evolve_group = hoomd.group.tag_list('tags',evolving_tag_list)

        ##group containing al beads that do not evolve
        #fixed_group = hoomd.group.tag_list('tags',fixed_tag_list)

        #defining walls at the boundaries to trap polymer fragments (see documentation at https://tinyurl.com/25xybm8c)
        walls = hoomd.md.wall.group()
        dh = 300
        top_zx = hoomd.md.wall.plane(origin=(0, 0, dh/2.), normal=(0, 0, -1), inside=True)
        bottom_zx = hoomd.md.wall.plane(origin=(0, 0, -dh/2.), normal=(0, 0, 1), inside=True)
        top_yz = hoomd.md.wall.plane(origin=(0, dh/2., 0), normal=(0, -1, 0), inside=True)
        bottom_yz = hoomd.md.wall.plane(origin=(0, -dh/2., 0), normal=(0, 1, 0), inside=True)
        top_xy = hoomd.md.wall.plane(origin=(dh/2., 0, 0), normal=(-1, 0, 0), inside=True)
        bottom_xy = hoomd.md.wall.plane(origin=(-dh/2., 0, 0), normal=(1, 0, 0), inside=True)
        walls.add(top_zx)
        walls.add(bottom_zx)
        walls.add(top_yz)
        walls.add(bottom_yz)
        walls.add(top_xy)
        walls.add(bottom_xy)
        wallLJ = hoomd.md.wall.lj(walls, r_cut = 3.);
        wallLJ.force_coeff.set('CA', sigma=1., epsilon=1.)

    #setting up the integrator
    hoomd.md.integrate.mode_standard(dt=@dt);
    integrator = hoomd.md.integrate.langevin(group = evolve_group, kT = 1, seed=@seed2*((1+icycle)**2)%20000)
    integrator.set_gamma('CA', gamma=1)
    

    ##log file to save some simulation data on the go
    #hoomd.analyze.log(filename=output_dir+'/log-output_'+fn_base_static+'.log',
    #              quantities=['potential_energy', 'temperature'],
    #              period=1000,
    #              overwrite=False);

    #Running simulation of DNA degradation with chopping at regular intervals. The simulation evolves only the beads included in evolve group.
    #The evolving beads become the pool of fragments that is defined by the removed bonds and angles.
    #A total of "ncut" bonds get chopped at random every cycle. "ncycles" is the total number of cycles.
    #The cuts are performed at the very beginning of a cycle and they accumulate as time passes.

    ##uncomment if equilibration is required
    #teq = @teq
    #hoomd.run(teq)
    
    if icycle == 0:
        #log file to save some simulation data on the go
        hoomd.analyze.log(filename=output_dir+'/log-output_'+fn_base_static+'.log',
                  quantities=['potential_energy', 'temperature'],
                  period=25000,
                  overwrite=False);
        
        #output file: all positions are saved by default, edit the argument of 'group=' to save whichever you want
        hoomd.dump.gsd(output_dir+'/degraded_'+fn_base_static+'.gsd', overwrite=False, truncate=False, period=@freq, group=evolve_group, phase=0)
          

    hoomd.run(nsteps)
    
    integrator.disable()
    
    #yy = system.take_snapshot(all=True)

    #saving the details free ends, kept bonds and kept angles at each cycle
    #also saving the list of chopped bonds and angles at the last cycle (at the very end)
    impt_data_file = output_dir+'/degradation_data_'+fn_base_static+'.npz'
    np.savez_compressed(impt_data_file,special_particles_record=special_particles_record, kept_bonds_record=kept_bonds_record, kept_angles_record=kept_angles_record, removed_bonds_record=removed_bonds_record,removed_angles_record=removed_angles_record,ncycles=ncycles,fin_bonds_record=fin_bonds_record,nsteps=nsteps,ncut=ncut) #,fixed_tag_list=fixed_tag_list)
    
    snapshot = system.take_snapshot(all=True)
    
#hoomd.dump.gsd(output_dir+'/test.gsd', overwrite=True,  period=100, group=all, phase=0)
