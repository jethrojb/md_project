# This module defines the class containing the positions and velocities of the particles

import numpy as np
import numba as nb

particle_spec = [('x', nb.float64), ('y', nb.float64), ('z', nb.float64),
                 ('vx', nb.float64), ('vy', nb.float64), ('vz', nb.float64), 
                 ('fx', nb.float64), ('fy', nb.float64), ('fz', nb.float64), 
                 ('dx', nb.float64), ('dy', nb.float64), ('dz', nb.float64), 
                 ('dr2', nb.float64), ('pe', nb.float64)]

sim_spec = [('T', nb.float64), ('rho', nb.float64), ('Nm', nb.int64), ('N', nb.int64),
            ('eq', nb.int64), ('pr', nb.int64), ('itrr', nb.int64), 
            ('rc', nb.float64), ('rc2', nb.float64), ('bl', nb.float64), 
            ('dt', nb.float64), ('k', nb.float64),
            ('length', nb.float64), ('output', nb.int64), ('utail', nb.float64), 
            ('ptail', nb.float64), ('seed', nb.float64), ('rdfmin', nb.float64), 
            ('rdfmax', nb.float64), ('rdfN', nb.int64), ('rdf', nb.int64),
            ('inputfile', nb.types.unicode_type), ('outputfile', nb.types.unicode_type),
            ('moviefile', nb.types.unicode_type), ('moviefreq', nb.int64),
            ('rescale_freq', nb.int64)]

@nb.experimental.jitclass(particle_spec)
class particle:
    def __init__(self):
        self.x=0.0              # x position
        self.y=0.0
        self.z=0.0
        self.vx=0.0             # x velocity
        self.vy=0.0
        self.vz=0.0
        self.fx=0.0             # x force
        self.fy=0.0
        self.fz=0.0
        self.dx=0.0             # x displacement for diffusion
        self.dy=0.0             # y displacement for diffusion
        self.dz=0.0             # z displacement for diffusion
        self.dr2=0.0            # MSD accumulator for diffusion
        self.pe=0.0

@nb.experimental.jitclass(sim_spec)
class simulation:
    def __init__(self):
        self.T=0.0              # temperature
        self.rho=0.0            # density
        self.Nm=0               # number of molecules
        self.N=0                # number of atoms
        self.eq=0               # equilibration steps
        self.pr=0               # production steps
        self.itrr=0             # interval for saving trajectories
        self.rc=0.0             # cutoff radius
        self.rc2=0.0            # square of cutoff radius
        self.bl=0.0             # bond length
        self.dt=0.0             # time step
        self.k=0.0
        self.length=0.0         # length of simulation blox
        self.output=0           # interval for output of instantaneous props
        self.utail=0.0          # tail correction to energy
        self.ptail=0.0          # tail correction to pressure
        self.seed=-1.0
        self.rdfmin=0.0         # minimum r value for rdf
        self.rdfmax=0.0         # maximum r value for rdf
        self.rdfN=0             # number of bins for rdf
        self.rdf=0              # frequency to accumulate the rdf
        self.inputfile=''       # name of input file
        self.outputfile=''      # name of output file
        self.moviefile=''
        self.moviefreq=0
        self.rescale_freq=0      # temperature rescale frequency during equilibration

class simprops:
    def __init__(self):
        self.ke=0.0             # kinetic energy
        self.pe=0.0             # potential energy
        self.pe2=0.0            # squared potential energy
        self.T=0.0              # temperature
        self.virial=0.0         # virial for pressure
        self.Nhist=0            # number of times accumulated
