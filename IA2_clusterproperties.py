"""
Script that calculates important cluster properties such as critical mass and radius.
Saves information to text as CGS units and human readable units.

Saved values are: r_Delta, M_Delta, rho_Delta, c_s, t_s.

Author: Angela Xue
Date: April 2026
"""

import numpy as np
from utils.units import CGS
from scipy.interpolate import RegularGridInterpolator

# Cosmological parameters (in CGS units).
Omega_r0, Omega_bm0, Omega_dm0, Omega_Lambda0, Omega_k0 = 0, 0.0441, 0.2139, 0.742, 0
h_0 = 0.72
H_0 = 100*h_0*1e5/(CGS.pc*1e6) # Converting to CGS units.

H2_a = lambda a : H_0**2 * (Omega_r0/a**4 + (Omega_bm0+Omega_dm0)/a**3 + Omega_Lambda0 + Omega_k0/a**2)
H2_z = lambda z : H2_a(1/(1+z))
rho_cr_z = lambda z: 3*H2_z(z)/(8*np.pi*CGS.G)

if __name__=="__main__":

    from scipy.optimize import root_scalar
    from utils.enzo import IA2, IA2Data, get_redshift
    import argparse, os
    
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("Delta", help="Delta times the critical density of the universe.", type=str)
    parser.add_argument("-snapshot", help="Which snapshot to calculate parameters of. Defaults to 15th.",required=False)
    parser.add_argument("-c", help="How small a subset to plot.",required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    Delta_list = [float(_Delta) for _Delta in args.Delta.split(',')]
    snapshot = int(args.snapshot) if args.snapshot else 15
    c = int(args.c) if args.c else 1

    # Prelim
    if 'E3A' in Path or 'E3A' in os.getcwd().split('/')[-1]:
        IA2.dim = 1024

    # Redshift constants
    z = get_redshift(snapshot)
    rho_cr = rho_cr_z(z)
    dL = IA2.dL*CGS.pc*1e3/(1+z) * c

    # Load data
    ia2 = IA2Data(Path)
    dens = np.array(ia2.get_val("Density",snapshot,c=c,units='cgs')+ia2.get_val("Dark_Matter_Density",snapshot,c=c,units='cgs'),dtype=np.float64)
    temp = np.array(ia2.get_val("Temperature",snapshot,c=c,units='cgs'),dtype=np.float64)

    # Calculate derived values
    def get_Rmesh():
        centre = ia2.get_centre(snapshot,c=c)
        xyz = tuple(np.arange(0,IA2.dim//c)-_centre for _centre in centre)
        XYZ = np.meshgrid(*xyz,indexing='ij')
        return np.sum([X**2 for X in XYZ],axis=0)**.5 * dL
    R = get_Rmesh()

    # Cluster properties.
    V_r = lambda r : 4*np.pi/3 * r**3
    r_arr = np.arange(0,R.max(),dL)
    M_arr = np.array([np.sum(dens[R<_r])*dL**3 for _r in r_arr])
    def M_r(r):
        try:
            return RegularGridInterpolator(r_arr[np.newaxis],M_arr)(r)
        except:
            return RegularGridInterpolator(r_arr[np.newaxis],M_arr)([r])[0]
    rhoavg_r = lambda r : M_r(r) / V_r(r)

    for Delta in Delta_list:
        r_Delta = root_scalar(lambda r : rhoavg_r(r)/rho_cr-Delta,x0=.1*(CGS.pc*1e6)).root
        M_Delta = M_r(r_Delta)
        rho_Delta = rhoavg_r(r_Delta)
        c_s = np.mean((IA2.gamma*CGS.kB*temp[R<r_Delta]/CGS.mp)**.5)
        t_s = 2*r_Delta/c_s

        with open(os.path.join(Path,f'cluster_{int(Delta)}.txt'),'w') as f:
            f.write(f"r_{int(Delta)} M_{int(Delta)} rho_{int(Delta)} c_s t_s\n")
            f.write(f"{r_Delta}, {M_Delta}, {rho_Delta}, {c_s}, {t_s}\n")
            f.write(f"{r_Delta/(CGS.pc*1e6)} Mpc, 10^{np.log10(M_Delta/CGS.Msun)} Msun, {rho_Delta/CGS.mp} m_p/cm^3, {c_s/1e5} km/s, {t_s/(CGS.yr*1e9)} Gyr\n")

    # np.savetxt(os.path.join(Path,f'cluster_{Delta}.txt'),[r_Delta,M_Delta,rho_Delta,c_s,t_s],header=f"r_{Delta} M_{Delta} rho_{Delta} c_s t_s")



