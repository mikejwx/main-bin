execfile('/home/xb899100/bin/projectFun.py')
execfile('/home/xb899100/bin/gma.py')
import scipy.io as io
data = io.readsav('/home/xb899100/Other/profilearmsonde_save_20010401-20060816.dat')
temp = data['tdry'][2,3:-1]
Q = data['qsat'][2,3:-1]*data['rh'][2,3:-1]/100.
q_kgkg = Q/1000.
temp_v = (temp+273.15)*(1 + 0.608*q_kgkg)
Z = data['alt'][2,3:-1]
p = data['pres'][3:-1]
myCAPE = getCAPE(T[0], T-273.15, Q, P, Z)

CAPE1 = 0
for level in range(len(Z)):
    if DT[level] > 0:
        dz = Z[level] - Z[level-1]
        CAPE1 += g*DT[level]*dz
def CAPEnew(T, Tp, P):
    lnp = np.log(P)
    Rd = 287.
    CAPE = 0
    for level in range(len(P)):
        dlnp = abs(lnp[level] - lnp[level-1])
        if Tp[level] > T[level]:
            CAPE += Rd*(Tp[level] - T[level])*dlnp
    
    return CAPE
def getCAPE2(T, Q, P, Z):
    theta_e = equivalent_potential_temperature(T-273.15, P, Q)
    theta_es = saturated_equivalent_potential_temperature(T-273.15, P, Q)
    CAPE = 0
    g = 9.81
    for level in range(len(Z)):
        if theta_e[0] > theta_es[level]:
            dz = Z[level] - Z[level-1]
            CAPE += g*(theta_e[0] - theta_es[level])*dz/theta_e[0]
    
    return CAPE