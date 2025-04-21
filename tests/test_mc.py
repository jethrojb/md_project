# This is to test that my code works

from src.diatomic_gas import particle, simulation, simprops
from src.initialize import initializeGrid, initializeVelocities, initializeFiles
from src.md import md
import numpy as np

sim = simulation()
sim.T = 8
sim.rho=0.0009
sim.Nm=100
sim.eq=10000
sim.pr=10000
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

initializeVelocities(sim, atom)

initializeFiles(sim, atom)

md(sim, atom)
