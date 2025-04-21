from src.rdf import rdf_finalize
import src.dhist as dh

def finalizefile(sim, atom, aprop, rdfcalls):
    # Variables
    pr=sim.pr
    N=sim.N
    
    # Calculate simple averages
    pe=aprop.pe/pr
    pe2=aprop.pe2/pr
    virial=aprop.virial/pr
    ke=aprop.ke/pr
    T=aprop.T/pr
    
    # Calculate heat capacity and pressure
    P = sim.rho*T + 1.0/3.0/sim.length**3.0*virial + sim.ptail
    cv = 3.0/2.0/(1 - 2.0/3.0*(pe2 - pe*pe)/N/(T*T))     # nve expression
   
    # Calculate the diffusivity from the MSD
    # This is zero for mc simulations
    Dmsd = 0.0
    for i in range(sim.Nm):
        Dmsd+=atom[i].dx*atom[i].dx + atom[i].dy*atom[i].dy + \
              atom[i].dz*atom[i].dz
    Dmsd=Dmsd/pr/N/6.0/sim.dt
    
    # Write the data to file
    fp=open(sim.outputfile, "a")
    
    # fp.write("\n    ***FINAL POSITIONS, XYZ Format***\n")
    # fp.write(str(N) + "\nYou can copy these coordinates to a file to " +
    #          "open in a viewer.\n")
    # for i in range(N):
    #     fp.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].x, \
    #                                                         atom[i].y, \
    #                                                         atom[i].z))

    # fp.write("\n         ***FINAL VELOCITIES***\n");
    # for i in range(N):
    #     fp.write("\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].vx, \
    #                                                     atom[i].vy, \
    #                                                     atom[i].vz))

    # fp.write("\n***Radial Distribution Function***\n\n");
    # rdf_finalize(sim, rdfh, rdfcalls)
    # dh.write(rdfh, fp)


    if sim.pr > 0:
        fp.write("\n***Simulation Averages***\n\n")
        fp.write("Temperature:            {:10.6f}\n".format(T))
        fp.write("Pressure:               {:10.6f}\n".format(P))
        fp.write("Heat Capacity:          {:10.6f}\n".format(cv))
        fp.write("Potential Energy:       {:10.6f}\n".format(pe/N + sim.utail))
        fp.write("Kinetic Energy:         {:10.6f}\n".format(ke/N))
        fp.write("Total Energy:           {:10.6f}\n".format((ke+pe) \
                                                                /N+sim.utail))
        fp.write("Diffusivity             {:10.6f}\n".format(Dmsd))
    else:
        fp.write("\nNo productions steps were specified, so simulation " +
                 "averages were not calculated.\n\n")
        
    fp.close()
