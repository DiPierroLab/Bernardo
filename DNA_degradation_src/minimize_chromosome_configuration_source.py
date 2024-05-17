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
# from hoomd import azplugins
import math

random.seed(@seed1)

hoomd.context.initialize("");

wallpos = @wallpos
snapshot = hoomd.data.make_snapshot(N=5500,
                                    box=hoomd.data.boxdim(Lx=wallpos, Ly=wallpos, Lz=wallpos),
                                    particle_types=['CA'],
                                    bond_types=['polymer'],
                                    angle_types=['backbone']);

input_traj = '@input_traj'
input_name = '@input_name'
loaded = md.load(input_traj)  # converted automatically to intended scaled positions, divided by 10 nm.
snap_index = @snap_index
pos = mb.load(loaded[int(snap_index)])
particle_pos = pos.xyz*10.0
xcom = np.mean(particle_pos,axis=0)
particle_pos = particle_pos - xcom
snapshot.particles.position[:] = particle_pos
n_particles = pos.n_particles
print(n_particles,particle_pos.max(axis=1))

snapshot.bonds.resize(n_particles-1)
i1 = 0
for i in range(n_particles-1):
    i2 = i1+1
    snapshot.bonds.group[i] = [i1,i2]
    i1 = i2

snapshot.angles.resize(n_particles-2)
i1 = 0
for i in range(n_particles-2):
    i2 = i1+1
    i3 = i2+1
    snapshot.angles.group[i] = [i1,i2,i3]
    i1 = i2

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
fene.bond_coeff.set('polymer',func=fene2, rmin=0, rmax=5.0, coeff = dict(kb=30, r0=1.5, rshift=1.5))
#fene.bond_coeff.set('polymer',func=fenewca, rmin=0.2, rmax=3.0, coeff = dict(kb=30, r0=1.5, rshift=1.5, epsilon=1.0, sigma=1.0))

angle = hoomd.md.angle.harmonic()
angle.angle_coeff.set('backbone',k=4,t0=np.pi)

lj = hoomd.md.pair.lj(r_cut=3.0, nlist=nl)
lj.pair_coeff.set('CA', 'CA', epsilon=1.0, sigma=1.0, r_cut=2**(1.0/6.0));
all = hoomd.group.all()

dt_for_fire = @dt_for_fire
maxirun = @maxirun
energy = []
fire=hoomd.md.integrate.mode_minimize_fire(dt=dt_for_fire, ftol=1e-2, Etol=1e-7) #ftol=1e-2, Etol=1e-7)
nve=hoomd.md.integrate.nve(group=all)
#irun = 0
#while not(fire.has_converged() or irun == maxirun):
#    print(fire.get_energy(), irun)
#    energy.append(fire.get_energy())
#    hoomd.run(1)
#    irun += 1

yy = system.take_snapshot()
xx = np.linalg.norm(np.diff(yy.particles.position,axis=0),axis=1)
while not(fire.has_converged()) and xx.min() < 1.0:
    hoomd.run(maxirun)
    yy = system.take_snapshot()
    xx = np.linalg.norm(np.diff(yy.particles.position,axis=0),axis=1)
    print(xx.min())

nve.disable()

dt_for_langvd = @dt_for_langvd
data_dir = '@data_dir'
hoomd.md.integrate.mode_standard(dt=dt_for_langvd);
integrator = hoomd.md.integrate.langevin(group = all, kT = 1, seed=@seed2)
integrator.set_gamma('CA', gamma=1)
hoomd.dump.gsd(data_dir+'/'+input_name+'_minimized.gsd',  period=1000, group=all, phase=0, overwrite=True)
hoomd.run(10000) 


