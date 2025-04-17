# This module contains code for initializing the simulation

# from src.diatomic_gas import 
import numpy as np
from src.diatomic_gas import particle
from numba.typed import List
from numba import njit


def initializeGrid(sim):
    """Accepts a simulation object and initializes grid positions."""
    sim.length = np.double(sim.N/sim.rho)**(1.0/3.0)
    sim.utail = (8/3)*np.pi*sim.rho*((1/3)*sim.rc**(-9.0) - sim.rc**(-3.0))
    sim.ptail = (16/3)*np.pi*sim.rho*sim.rho*((2/3)*sim.rc**(-9.0) - sim.rc**(-3.0))

    atom=List()

    atomindex=0
    case=0

    nlin = np.uint((sim.Nm/4)**(1.0/3.0))
    if nlin**3 < sim.Nm/4:
        nlin = nlin + 1

    a = sim.length/nlin

    for z in range(nlin):
        for y in range(nlin):
            for x in range(nlin):
                for i in range(4):
                    if atomindex == sim.Nm:
                        return atom
                    
                    atom.append(particle())
                    if case == 0:
                        atom[atomindex].x = 0.0+x*a
                        atom[atomindex].y = 0.0+y*a
                        atom[atomindex].z = 0.0+z*a
                        case=1
                    elif case == 1:
                        atom[atomindex].x = 0.0+x*a
                        atom[atomindex].y = 0.5*a+y*a
                        atom[atomindex].z = 0.5*a+z*a
                        case=2
                    elif case == 2:
                        atom[atomindex].x = 0.5*a+x*a
                        atom[atomindex].y = 0.0+y*a
                        atom[atomindex].z = 0.5*a+z*a
                        case=3
                    else:
                        atom[atomindex].x = 0.5*a+x*a
                        atom[atomindex].y = 0.5*a+y*a
                        atom[atomindex].z = 0.0+z*a
                    
                    atomindex += 1

    return atom
                       





    