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

hoomd.context.initialize("");
minimized_str_file = '@minimized_str_file'
f = gsd.fl.open(minimized_str_file,'rb')

snapshot = hoomd.data.gsd_snapshot(minimized_str_file,frame=f.nframes-1)
n_particles = snapshot.particles.N

random.seed(@seed1)

nc = @nkeep # number of bonds to keep after chopping
all_bond_list = np.arange(n_particles-1)
bond_list = np.sort(np.random.choice(range(n_particles-1), nc, replace=False))
removed_bond_list = np.setdiff1d(all_bond_list,bond_list)
nr = len(removed_bond_list)
snapshot.bonds.resize(nc)

all_bonds = []
i1 = 0
for i in range(n_particles-1):
    i2 = i1+1
    all_bonds.append([i1,i2])
    i1 = i2

for i in range(nc):
    snapshot.bonds.group[i] = all_bonds[bond_list[i]]

special_particles = []
for i in range(len(removed_bond_list)):
    m = removed_bond_list[i]
#     print(all_bonds[m])
    special_particles.append(all_bonds[m][0])
    special_particles.append(all_bonds[m][1])    

special_particles = np.unique(special_particles)
n_sp = len(special_particles)

c = @chr  # chromosome identity
r = '@region'  # region identity
filenumber = @filenumber
str_no = @str
run_no = @run
output_dir = '@output_dir'
fn_base = 'C_'+str(c)+'_region_'+str(r)+'_nkeep_'+str(nc)+'_str_'+str(str_no)+'_n_'+str(filenumber)+'_sim_'+str(run_no)

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

system = hoomd.init.read_snapshot(snapshot);

nl = hoomd.md.nlist.cell();
nl.reset_exclusions(exclusions = []);

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

fene = hoomd.md.bond.table(width=1000)
fene.bond_coeff.set('polymer',func=fene2, rmin=0., rmax=5.0, coeff = dict(kb=30, r0=1.5, rshift=1.5))

angle = hoomd.md.angle.harmonic()
angle.angle_coeff.set('backbone',k=4,t0=np.pi)

lj = hoomd.md.pair.lj(r_cut=3.0, nlist=nl)
lj.pair_coeff.set('CA', 'CA', epsilon=1.0, sigma=1.0, r_cut=2**(1.0/6.0));

all = hoomd.group.all()

n = @nevolve
all_list = np.arange(n_particles)
evolving_tag_list = np.sort(np.random.choice(range(n_particles), n, replace=False))
fixed_tag_list = np.setdiff1d(all_list,evolving_tag_list)
# np.sort(fixed_tag_list),len(fixed_tag_list)


impt_data_file = output_dir+'/chopping_and_ligation_data_'+fn_base+'.npz'
np.savez_compressed(impt_data_file,ligating_ends=special_particles,removed_bonds=removed_bond_list,removed_angles=removed_angles,fixed_tag_list=fixed_tag_list)

evolve_group = hoomd.group.tag_list('tags',evolving_tag_list)

fixed_group = hoomd.group.tag_list('tags',fixed_tag_list)

walls = hoomd.md.wall.group()
dh = @wallpos
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

hoomd.md.integrate.mode_standard(dt=@dt);

integrator = hoomd.md.integrate.langevin(group = evolve_group, kT = 1, seed=@seed2)
integrator.set_gamma('CA', gamma=1)

hoomd.analyze.log(filename=output_dir+'/log-output_'+fn_base+'.log',
                  quantities=['potential_energy', 'temperature'],
                  period=5000,
                  overwrite=True);

teq = @teq
nsteps = @nsteps

hoomd.run(teq)

hoomd.dump.gsd(output_dir+'/chromosome_crosslinked_fixed_chopped_'+fn_base+'.gsd', overwrite=True,  period=@freq, group=all, phase=0)
#hoomd.dump.gsd('chromosome_evolving_fixed_chopped.gsd', overwrite=True,  period=5000, group=evolve_group, phase=0)
#hoomd.dump.gsd('chromosome_fixed_chopped.gsd', overwrite=True,  period=5000, group=fixed_group, phase=0)

hoomd.run_upto(nsteps)

