from src.forces import forces
from src.kinetic import ke_and_T
from src.verlet import verlet1, verlet2
from src.diatomic_gas import simprops
from src.scale_velocities import scalevelocities
import numpy as np
import src.dhist as dh
from src.rdf import rdf_accumulate
from src.finalize_file import finalizefile

def md(sim, atom):
    rescale_freq = sim.rescale_freq

    instant_prop = simprops()
    avg_prop = simprops()

    instant_prop.pe, instant_prop.virial = forces(sim, atom)
    instant_prop.ke, instant_prop.T = ke_and_T(atom)
    P = sim.rho*instant_prop.T + 1.0/3.0/sim.length**3.0*instant_prop.virial + sim.ptail
    fp = open(sim.outputfile, "a")
    fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}    {:13.6f}    " \
             "{:13.6f}    {:13.6f}    {:13.6f}\n" \
             .format(0, instant_prop.T, instant_prop.T, P, P, instant_prop.ke/sim.N, \
                     instant_prop.pe/sim.N + sim.utail, \
                     (instant_prop.ke + instant_prop.pe)/sim.N + sim.utail))
    fp.close()


    for i in range(1, np.int64(sim.eq+1)):
        verlet1(sim, atom)
        instant_prop.pe, instant_prop.virial = forces(sim, atom)
        verlet2(sim, atom)
        instant_prop.ke, instant_prop.T = ke_and_T(atom)

        avg_prop.pe = avg_prop.pe + instant_prop.pe
        avg_prop.ke = avg_prop.ke + instant_prop.ke
        avg_prop.T = avg_prop.T + instant_prop.T
        avg_prop.virial = avg_prop.virial + instant_prop.virial

        if i % sim.output == 0:
            P=sim.rho*instant_prop.T + 1.0/3.0/sim.length**3.0*instant_prop.virial + \
              sim.ptail
            Pave=sim.rho*avg_prop.T/i + 1.0/3.0/sim.length**3.0*avg_prop.virial/i + \
              sim.ptail 
            fp=open(sim.outputfile, "a")
            fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}    {:13.6f}    " \
                     "{:13.6f}    {:13.6f}    {:13.6f}\n" \
                     .format(i, instant_prop.T, avg_prop.T/i, P, Pave, instant_prop.ke/sim.N, \
                             instant_prop.pe/sim.N + sim.utail, \
                             (instant_prop.ke + instant_prop.pe)/sim.N + sim.utail))
            fp.close()
            print("Equilibration Step " + str(i) + "\n")

        if i % rescale_freq == 0:
            scalevelocities(sim, atom, avg_prop.T/i)

    avg_prop.pe = 0.0
    avg_prop.ke = 0.0
    avg_prop.T = 0.0
    avg_prop.virial = 0.0

    for i in range(sim.N):
        atom[i].dx = 0.0
        atom[i].dy = 0.0
        atom[i].dz = 0.0
    
    Nrdfcalls = 0
    # rdfh = dh.hist(sim.rdfmin, sim.rdfmax, sim.rdfN)

    for i in range(1, np.int64(sim.pr+1)):
        verlet1(sim, atom)
        instant_prop.pe, instant_prop.virial = forces(sim, atom)
        verlet2(sim, atom)
        instant_prop.ke, instant_prop.T = ke_and_T(atom)

        avg_prop.pe += instant_prop.pe
        avg_prop.ke += instant_prop.ke
        avg_prop.T += instant_prop.T
        avg_prop.virial += instant_prop.virial
        avg_prop.pe2 += instant_prop.pe*instant_prop.pe

        if i % sim.output == 0:
            P=sim.rho*instant_prop.T + 1.0/3.0/sim.length**3.0*instant_prop.virial + \
              sim.ptail
            Pave=sim.rho*avg_prop.T/i + 1.0/3.0/sim.length**3.0*avg_prop.virial/i + \
              sim.ptail 
            fp=open(sim.outputfile, "a")
            fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}    {:13.6f}    " \
                     "{:13.6f}    {:13.6f}    {:13.6f}\n" \
                     .format(i, instant_prop.T, avg_prop.T/i, P, Pave, instant_prop.ke/sim.N, \
                             instant_prop.pe/sim.N + sim.utail, \
                             (instant_prop.ke + instant_prop.pe)/sim.N + sim.utail))
            fp.close()
            print("Production Step " + str(i) + "\n")

            # if i % sim.rdf == 0:
                # Nrdfcalls += 1
                # rdf_accumulate(sim, atom, rdfh)
        if sim.moviefreq != 0:
          if i % sim.moviefreq == 0:
              fm = open(sim.moviefile, 'a')
              fm.write(str(sim.N)+'\n\n')
              for i in range(sim.N):
                      fm.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].x, \
                                                                          atom[i].y, \
                                                                          atom[i].z))
              fm.close()

    finalizefile(sim, atom, avg_prop, Nrdfcalls)

