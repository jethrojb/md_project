# This is to test that my code works

from src.diatomic_gas import particle, simulation, simprops
import numpy as np

sim = simulation()
sim.T = 0.85
sim.rho=0.9
sim.Nm=50
sim.eq=50
sim.pr=100
sim.itrr=1
sim.rc=3.0
sim.dt=0.01
sim.rdfmin=0.8
sim.rdfmax=4.0
sim.sdfN=100
sim.rdf=1

