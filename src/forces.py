# import numpy as np
# from numba import njit

# @njit
# def forces(sim, atom):
#     hL = np.float64(sim.length*0.5)
#     N = np.int64(sim.N)
#     Nm = np.int64(sim.Nm)

#     for i in range(sim.N):
#         atom[i].fx = 0.0
#         atom[i].fy = 0.0
#         atom[i].fz = 0.0
    
#     virial = 0.0
#     pe = 0.0

#     for i in range(N - 1):
#         for j in range(i+1, N):
#             # skip atoms that are bonded
#             if j == i + sim.Nm:
#                 continue

#             dx = atom[i].x - atom[j].x
#             dy = atom[i].y - atom[j].y
#             dz = atom[i].z - atom[j].z

#             if np.abs(dx) > hL:
#                 if dx < 0.0:
#                     dx = dx + sim.length
#                 else:
#                     dx = dx - sim.length
#             if np.abs(dy) > hL:
#                 if dy < 0.0:
#                     dy = dy + sim.length
#                 else:
#                     dy = dy - sim.length
#             if np.abs(dz) > hL:
#                 if dz < 0.0:
#                     dz = dz + sim.length
#                 else:
#                     dz = dz - sim.length

#             dr2 = dx*dx + dy*dy + dz*dz

#             if dr2 < sim.rc2:
#                 d2 = 1.0/dr2
#                 d4 = d2*d2
#                 d8 = d4*d4
#                 d14 = d8*d4*d2

#                 fr = 48.0*(d14 - 0.5*d8)

#             # components of forces
#                 atom[i].fx=atom[i].fx+fr*dx
#                 atom[i].fy=atom[i].fy+fr*dy
#                 atom[i].fz=atom[i].fz+fr*dz
#                 atom[j].fx=atom[j].fx-fr*dx
#                 atom[j].fy=atom[j].fy-fr*dy
#                 atom[j].fz=atom[j].fz-fr*dz

#                 virial += dr2*fr
#                 pe += 4.0*(d14 - d8)*dr2

#     # Now we need to deal with bonded atoms
#     # for i in range(Nm):
#     #     # Ensure i and i + Nm are valid indices
#     #     # This loop structure assumes atoms 0 to Nm-1 are first atoms
#     #     # and atoms Nm to 2*Nm-1 are second atoms.

#     #     dx = atom[i].x - atom[i + Nm].x
#     #     dy = atom[i].y - atom[i + Nm].y
#     #     dz = atom[i].z - atom[i + Nm].z

#     #     # Apply periodic boundary conditions to the bond vector too
#     #     # This is important for calculating the correct distance across boundaries
#     #     if np.abs(dx) > hL:
#     #         if dx < 0.0:
#     #             dx = dx + sim.length
#     #         else:
#     #             dx = dx - sim.length
#     #     if np.abs(dy) > hL:
#     #         if dy < 0.0:
#     #             dy = dy + sim.length
#     #         else:
#     #             dy = dy - sim.length
#     #     if np.abs(dz) > hL:
#     #         if dz < 0.0:
#     #             dz = dz + sim.length
#     #         else:
#     #             dz = dz - sim.length


#     #     dr2 = dx*dx + dy*dy + dz*dz
#     #     dr = np.sqrt(dr2)

#     #     # Avoid division by zero if dr is somehow zero (shouldn't happen for bonded atoms)
#     #     if dr > 1e-6: # Use a small tolerance
#     #         # Calculate the magnitude of the force, divided by the distance dr
#     #         # The force magnitude is -sim.k * (dr - sim.bl)
#     #         # The force vector is proportional to (dx, dy, dz) / dr
#     #         force_over_r = -sim.k * (dr - sim.bl) / dr

#     #         # Components of forces
#     #         # Force on atom i points away from atom i+Nm if dr < sim.bl
#     #         # Force on atom i points towards atom i+Nm if dr > sim.bl
#     #         # The vector from i+Nm to i is (dx, dy, dz)
#     #         # So force on i is proportional to (dx, dy, dz) * force_over_r
#     #         atom[i].fx = atom[i].fx + force_over_r * dx
#     #         atom[i].fy = atom[i].fy + force_over_r * dy
#     #         atom[i].fz = atom[i].fz + force_over_r * dz

#     #         # Force on atom i + Nm is equal and opposite
#     #         atom[i + Nm].fx = atom[i + Nm].fx - force_over_r * dx
#     #         atom[i + Nm].fy = atom[i + Nm].fy - force_over_r * dy
#     #         atom[i + Nm].fz = atom[i + Nm].fz - force_over_r * dz

#     #         # Virial contribution from this bond
#     #         virial += dr2 * force_over_r # This is r^2 * (F/r) = r * F
#     #                                      # So this adds -sim.k * (dr - sim.bl) * dr2 / dr = -sim.k * (dr - sim.bl) * dr
#     #                                      # This is the correct virial contribution for a harmonic bond

#     #         # Potential energy contribution from this bond
#     #         pe += 0.5 * sim.k * (dr - sim.bl) * (dr - sim.bl)
#     # # Now we need to deal with bonded atoms
#     # for i in range(Nm):
#     #     dx = atom[i].x - atom[i + Nm].x
#     #     dy = atom[i].y - atom[i + Nm].y
#     #     dz = atom[i].z - atom[i + Nm].z

#     #     dr2 = dx*dx + dy*dy + dz*dz
#     #     dr = np.sqrt(dr2)

#     #     fr = -sim.k*(dr - sim.bl)

#     #     atom[i].fx = atom[i].fx + fr*dx
#     #     atom[i].fy = atom[i].fy + fr*dy
#     #     atom[i].fz = atom[i].fz + fr*dz
#     #     atom[i + Nm].fx = atom[i + Nm].fx - fr*dx
#     #     atom[i + Nm].fy = atom[i + Nm].fy - fr*dy
#     #     atom[i + Nm].fz = atom[i + Nm].fz - fr*dz

#     #     virial += dr2*fr
#     #     pe += 0.5*sim.k*(dr - sim.bl)*(dr - sim.bl)

#     return (np.float64(pe), np.float64(virial))

import numpy as np
from numba import njit

# This function is passed a simulation object and a list of site object.
# It returns the potential energy of the system and also assigns the 
# forces on each site. 
@njit
def forces(sim,atom):
    # Variables
    hL=np.float64(sim.length*0.5)   # half the box length
    N=np.int64(sim.N)
    
    # Zero out the force accumulators for each particle
    for i in range(sim.N):
        atom[i].fx=0.0
        atom[i].fy=0.0
        atom[i].fz=0.0
    
    # Zero out the system accumulators
    virial=0.0  # virial portion of pressure
    pe=0.0      # potential energy
       
    # Calculate the forces by looping over all pairs of sites
    for i in range(N-1):
        for j in range(i+1, N):
            # Calculate the distance between sites i and j
            dx=atom[i].x-atom[j].x
            dy=atom[i].y-atom[j].y
            dz=atom[i].z-atom[j].z
            
            # Minimum image convention
            if np.abs(dx)>hL:
                if dx < 0.0: dx=dx+sim.length
                else: dx=dx-sim.length
            if np.abs(dy)>hL:
                if dy < 0.0: dy=dy+sim.length
                else: dy=dy-sim.length
            if np.abs(dz)>hL:
                if dz < 0.0: dz=dz+sim.length
                else: dz=dz-sim.length
            
            dr2=dx*dx+dy*dy+dz*dz
            
            # Calculate the energy and force for the pair
            if dr2 < sim.rc2: # apply cutoff
                d2=1.0/dr2
                d4=d2*d2
                d8=d4*d4
                d14=d8*d4*d2
                fr=48.0*(d14-0.5*d8)
                
                
                # components of forces
                atom[i].fx=atom[i].fx+fr*dx
                atom[i].fy=atom[i].fy+fr*dy
                atom[i].fz=atom[i].fz+fr*dz
                atom[j].fx=atom[j].fx-fr*dx
                atom[j].fy=atom[j].fy-fr*dy
                atom[j].fz=atom[j].fz-fr*dz   
                
                # virial and potential energy
                virial=virial+dr2*fr
                pe=pe+4.0*(d14-d8)*dr2

    return(np.float64(pe), np.float64(virial))

