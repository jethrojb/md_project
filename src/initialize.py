# This module contains code for initializing the simulation

# from src.diatomic_gas import 
import numpy as np
from src.diatomic_gas import particle, simulation
from numba.typed import List
from numba import njit


def initializeGrid(sim):
    """Accepts a simulation object and initializes grid positions."""
    sim.length = np.double(sim.Nm/sim.rho)**(1.0/3.0)
    sim.utail = (8/3)*np.pi*sim.rho*((1/3)*sim.rc**(-9.0) - sim.rc**(-3.0))
    sim.ptail = (16/3)*np.pi*sim.rho*sim.rho*((2/3)*sim.rc**(-9.0) - sim.rc**(-3.0))
    sim.N = 2*sim.Nm

    atom=List()
    for i in range(2*sim.Nm):
        atom.append(particle())


    moleculeindex=0
    case=0

    nlin = np.uint((sim.Nm/4)**(1.0/3.0))

    if nlin**3 < sim.Nm/4:
        nlin = nlin + 1

    a = sim.length/nlin

    for z in range(nlin):
        for y in range(nlin):
            for x in range(nlin):
                for i in range(4):
                    if moleculeindex == sim.Nm:
                        return atom
                    if case == 0:
                        atom[moleculeindex].x = 0.0+x*a
                        atom[moleculeindex].y = 0.0+y*a
                        atom[moleculeindex].z = 0.0+z*a
                        atom[moleculeindex + sim.Nm].x = 0.0+x*a + sim.bl
                        atom[moleculeindex + sim.Nm].y = 0.0+y*a
                        atom[moleculeindex + sim.Nm].z = 0.0+z*a

                        case=1
                    elif case == 1:
                        atom[moleculeindex].x = 0.0+x*a
                        atom[moleculeindex].y = 0.5*a+y*a
                        atom[moleculeindex].z = 0.5*a+z*a
                        atom[moleculeindex + sim.Nm].x = 0.0+x*a + sim.bl
                        atom[moleculeindex + sim.Nm].y = 0.5*a+y*a
                        atom[moleculeindex + sim.Nm].z = 0.5*a+z*a

                        case=2
                    elif case == 2:
                        atom[moleculeindex].x = 0.5*a+x*a
                        atom[moleculeindex].y = 0.0+y*a
                        atom[moleculeindex].z = 0.5*a+z*a
                        atom[moleculeindex + sim.Nm].x = 0.5*a+x*a + sim.bl
                        atom[moleculeindex + sim.Nm].y = 0.0+y*a
                        atom[moleculeindex + sim.Nm].z = 0.5*a+z*a

                        case=3
                    else:
                        atom[moleculeindex].x = 0.5*a+x*a
                        atom[moleculeindex].y = 0.5*a+y*a
                        atom[moleculeindex].z = 0.0+z*a
                        atom[moleculeindex + sim.Nm].x = 0.5*a+x*a + sim.bl
                        atom[moleculeindex + sim.Nm].y = 0.5*a+y*a
                        atom[moleculeindex + sim.Nm].z = 0.0+z*a

                    
                    moleculeindex += 1

    return atom
                       





    