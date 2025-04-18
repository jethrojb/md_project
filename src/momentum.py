def checkMomentum(atom):
    N = len(atom)
    
    vcumx = 0.0
    vcumy = 0.0
    vcumz = 0.0

    for i in range(N):
        vcumx = vcumx + atom[i].vx
        vcumy = vcumy + atom[i].vy
        vcumz = vcumz + atom[i].vz
    if(vcumx+vcumy+vcumz)< 10.0**-10.0: return(0)
    else: return(1)

def momentumCorrection(atom):
    N = len(atom)

    vcumx = 0.0
    vcumy = 0.0
    vcumz = 0.0

    for i in range(N):
        vcumx = vcumx + atom[i].vx
        vcumy = vcumy + atom[i].vy
        vcumz = vcumz + atom[i].vz

    vcumx = vcumx/N
    vcumy = vcumy/N
    vcumz /= N

    for i in range(N):
        atom[i].vx = atom[i].vx - vcumx
        atom[i].vy = atom[i].vy - vcumy
        atom[i].vz = atom[i].vz - vcumz

    return checkMomentum(atom)
