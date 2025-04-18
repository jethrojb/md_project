# This is to test that my code works

from src.diatomic_gas import particle, simulation, simprops
from src.initialize import initializeGrid, initializeVelocities, initializeFiles
from src.md import md
import numpy as np

sim = simulation()
sim.T = 0.85
sim.rho=0.9
sim.Nm=27
sim.eq=50
sim.pr=100
sim.itrr=1
sim.rc=3.0
sim.bl=0.329
sim.dt=0.01
sim.rdfmin=0.8
sim.rdfmax=4.0
sim.rdfN=100
sim.rdf=1
sim.outputfile='test_output.txt'
sim.output=1
sim.k = 488331.78

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
