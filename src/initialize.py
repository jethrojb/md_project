# This module contains code for initializing the simulation

# from src.diatomic_gas import 
import random, sys
import numpy as np
from src.diatomic_gas import particle, simulation
from src.kinetic import temperature
from src.scale_velocities import scalevelocities
from src.momentum import checkMomentum, momentumCorrection
from numba.typed import List
from numba import njit


def initializeGrid(sim):
    """Accepts a simulation object and initializes grid positions."""
    sim.length = np.double(sim.Nm/sim.rho)**(1.0/3.0)
    sim.utail = (8/3)*np.pi*sim.rho*((1/3)*sim.rc**(-9.0) - sim.rc**(-3.0))
    sim.ptail = (16/3)*np.pi*sim.rho*sim.rho*((2/3)*sim.rc**(-9.0) - sim.rc**(-3.0))
    sim.N = 2*sim.Nm
    sim.rc2 = sim.rc*sim.rc

    atom=List()
    for i in range(2*sim.Nm):
        atom.append(particle())


    moleculeindex=0
    case=0

    nlin = np.uint((sim.Nm/4.0)**(1.0/3.0))

    if nlin**3 < sim.Nm/4.0:
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

                        case=0

                    moleculeindex += 1

    return atom

def initializeVelocities(sim, atom):
    for i in range(sim.Nm):
        vx = random.uniform(-1,1)
        vy = random.uniform(-1,1)
        vz = random.uniform(-1,1)
        atom[i].vx = vx
        atom[i + sim.Nm].vx = vx
        atom[i].vy = vy
        atom[i + sim.Nm].vy = vy
        atom[i].vz = vz
        atom[i + sim.Nm].vz = vz

    momentumflag = momentumCorrection(atom)
    if momentumflag:
        sys.exit("The linear momentum could not be zeroed.")

    T = temperature(atom)

    scalevelocities(sim, atom, T)

def initializeFiles(sim, atom):
    # Initilaize the .trr file.
    fi=open(sim.moviefile,"w")
    fi.close()
    
    # Write the .xyz file needed for loading the .trr file.
    fi=open(sim.moviefile,"w")
    fi.write(str(sim.N)+"\nLoad this file in VMD\n")
    for i in range(sim.N):
        fi.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].x, \
                                                            atom[i].y, \
                                                            atom[i].z))
    fi.close()

    file = open(sim.outputfile, 'w')
    file.write('MD simulation of ' + str(sim.Nm) + 'LJ Particles at T*={:.4f}'.format(sim.T) +
               ' and ' + 'rho*={:.4f}\n\n'.format(sim.rho))
    
    file.write('Output File:        ' + sim.outputfile + '\n\n')
    file.write("\n    ***Input Parameters***\n");
    file.write("N           " + str(sim.Nm) + "\n")
    file.write("temp        " + str(sim.T) + "\n")
    file.write("rho         " + str(sim.rho) + "\n")
    file.write("esteps      " + str(sim.eq) + "\n")
    file.write("psteps      " + str(sim.pr) + "\n")
    file.write("rcut        " + str(sim.rc) + "\n")
    file.write("dt          " + str(sim.dt) + "\n")

    # file.write("\n    ***INITIAL POSITIONS, XYZ Format***\n")
    # file.write(str(sim.Nm) + "\nYou can copy these coordinates to a file to " +
    #          "open in a viewer.\n")
    # for i in range(sim.Nm):
    #     file.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].x, \
    #                                                         atom[i].y, \
    #                                                         atom[i].z))
    
    file.write("\n         ***INITIAL VELOCITIES***\n")
    for i in range(sim.Nm):
        file.write("\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].vx, \
                                                               atom[i].vy, \
                                                               atom[i].vz))
        
    file.write("\n\nIteration                T              T Ave.       " +
                 "       P              P Ave.            KE               " +
                 "PE               TE\n\n")
    
    file.close()
