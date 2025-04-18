from numba import njit
# This function is passed a simulation object and a list of site objects.
# It is the first needed to use the velocity verlet algorithm. It uses
# the data at time step t to update the positions to the next time step and
# the velocities to the next half time step.
@njit
def verlet1(sim, atom):
    
    for i in range(sim.N):
        # Update the positions to a full time step.
        dx=sim.dt*atom[i].vx+sim.dt*sim.dt*atom[i].fx/2.0
        dy=sim.dt*atom[i].vy+sim.dt*sim.dt*atom[i].fy/2.0
        dz=sim.dt*atom[i].vz+sim.dt*sim.dt*atom[i].fz/2.0
        atom[i].x=atom[i].x+dx
        atom[i].y=atom[i].y+dy
        atom[i].z=atom[i].z+dz
        
        # Update the displacement accumulators for diffusivity
        atom[i].dx=atom[i].dx+dx
        atom[i].dy=atom[i].dy+dy
        atom[i].dz=atom[i].dz+dz
        
        # Apply periodic boundary conditions
        if atom[i].x < 0.0:          atom[i].x=atom[i].x+sim.length
        elif atom[i].x > sim.length: atom[i].x=atom[i].x-sim.length
        if atom[i].y < 0.0:          atom[i].y=atom[i].y+sim.length
        elif atom[i].y > sim.length: atom[i].y=atom[i].y-sim.length
        if atom[i].z < 0.0:          atom[i].z=atom[i].z+sim.length
        elif atom[i].z > sim.length: atom[i].z=atom[i].z-sim.length
        
        # Update the velocities to half a time step
        atom[i].vx=atom[i].vx+sim.dt*atom[i].fx/2.0
        atom[i].vy=atom[i].vy+sim.dt*atom[i].fy/2.0
        atom[i].vz=atom[i].vz+sim.dt*atom[i].fz/2.0
        
# This function is passed a simulation object and a list of site objects.
# It is the second function needed to use the velocity verlet       
# algorithm.  It updates the velocites from the half time step to the 
# full time step.
@njit
def verlet2(sim, atom):
    for i in range(sim.N):
        atom[i].vx=atom[i].vx+sim.dt*atom[i].fx/2.0
        atom[i].vy=atom[i].vy+sim.dt*atom[i].fy/2.0
        atom[i].vz=atom[i].vz+sim.dt*atom[i].fz/2.0
