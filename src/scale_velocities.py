import numpy as np
def scalevelocities(sim, atom, temp):
    scale = np.sqrt(sim.T/temp)
    for i in range(sim.N):
        atom[i].vx *= scale
        atom[i].vy *= scale
        atom[i].vz *= scale
