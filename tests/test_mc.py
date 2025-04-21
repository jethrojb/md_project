# This is to test that my code works

from src.diatomic_gas import particle, simulation, simprops
from src.initialize import initializeGrid, initializeVelocities, initializeFiles
from src.md import md
import numpy as np

sim = simulation()
sim.T = 8
sim.rho=0.0009
sim.Nm=100
sim.eq=1000
sim.pr=1000
sim.itrr=1
sim.rc=3.0
sim.bl=0.329
sim.dt=0.01
sim.rdfmin=0.8
sim.rdfmax=4.0
sim.rdfN=100
sim.rdf=1
sim.outputfile='test_output.txt'
sim.output=100
sim.k = 0.8
sim.moviefile='movie.xyz'
sim.moviefreq=10
sim.rescale_freq=10

atom = initializeGrid(sim)

# print(atom)

# Assuming 'atom' is your numba.typed.List
# atom_array = np.array([p for p in atom])
# for i in range(sim.Nm):
#     print('x')
#     print(atom_array[i].x, atom_array[i + sim.Nm].x)
#     print('y')
#     print(atom_array[i].y, atom_array[i + sim.Nm].y)


# Test initializing velocities

initializeVelocities(sim, atom)

# Test initializing output file

initializeFiles(sim, atom)

md(sim, atom)
