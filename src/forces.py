import numpy as np
from numba import njit

@njit
def forces(sim, atom):
    hL = np.float64(sim.length*0.5)
    N = np.int64(sim.N)
    Nm = np.int64(sim.Nm)

    for i in range(sim.N):
        atom[i].fx = 0.0
        atom[i].fy = 0.0
        atom[i].fz = 0.0
    
    virial = 0.0
    pe = 0.0

    for i in range(N - 1):
        for j in range(i+1, N):
            # skip atoms that are bonded
            if j == i + sim.Nm:
                continue

            dx = atom[i].x - atom[j].x
            dy = atom[i].y - atom[j].y
            dz = atom[i].z - atom[j].z

            if np.abs(dx) > hL:
                if dx < 0.0:
                    dx = dx + sim.length
                else:
                    dx = dx - sim.length
            if np.abs(dy) > hL:
                if dy < 0.0:
                    dy = dy + sim.length
                else:
                    dy = dy - sim.length
            if np.abs(dz) > hL:
                if dz < 0.0:
                    dz = dz + sim.length
                else:
                    dz = dz - sim.length

            dr2 = dx*dx + dy*dy + dz*dz

            if dr2 < sim.rc2:
                d2 = 1.0/dr2
                d4 = d2*d2
                d8 = d4*d4
                d14 = d8*d4*d2

                fr = 48.0*(d14 - 0.5*d8)

            # components of forces
                atom[i].fx=atom[i].fx+fr*dx
                atom[i].fy=atom[i].fy+fr*dy
                atom[i].fz=atom[i].fz+fr*dz
                atom[j].fx=atom[j].fx-fr*dx
                atom[j].fy=atom[j].fy-fr*dy
                atom[j].fz=atom[j].fz-fr*dz

                virial += dr2*fr
                pe += 4.0*(d14 - d8)*dr2

    # Now we need to deal with bonded atoms
    for i in range(Nm):
        dx = atom[i].x - atom[i + Nm].x
        dy = atom[i].y - atom[i + Nm].y
        dz = atom[i].z - atom[i + Nm].z

        dr2 = dx*dx + dy*dy + dz*dz
        dr = np.sqrt(dr2)

        fr = -sim.k*(dr - sim.bl)

        atom[i].fx = atom[i].fx + fr*dx
        atom[i].fy = atom[i].fy + fr*dy
        atom[i].fz = atom[i].fz + fr*dz
        atom[i + Nm].fx = atom[i + Nm].fx - fr*dx
        atom[i + Nm].fy = atom[i + Nm].fy - fr*dy
        atom[i + Nm].fz = atom[i + Nm].fz - fr*dz

        virial += dr2*fr
        pe += 0.5*sim.k*(dr - sim.bl)*(dr - sim.bl)

    return (np.float64(pe), np.float64(virial))


