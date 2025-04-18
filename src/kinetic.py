def kinetic_energy(atom):
    N = len(atom)
    
    ke = 0.0
    
    for i in range(N):
        v2 = atom[i].vx*atom[i].vx + atom[i].vy*atom[i].vy + atom[i].vz*atom[i].vz
        ke = ke + 0.5*v2
    
    return ke
    
def temperature(atom):
    N = len(atom)

    T = (2.0/3.0/N)*kinetic_energy(atom)
    return T
