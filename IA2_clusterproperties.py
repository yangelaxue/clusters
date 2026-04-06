"""
Script that calculates important cluster properties such as critical mass and radius.

Author: Angela Xue
Date: April 2026
"""

import numpy as np
from utils.units import CGS
from scipy.interpolate import RegularGridInterpolator

Omega_r0, Omega_bm0, Omega_dm0, Omega_Lambda0, Omega_k0 = 0, 0.0441, 0.2139, 0.742, 0
h_0 = 0.72
H_0 = 100*h_0*1e5/(CGS.pc*1e6) # Converting to CGS units.

# Hubble factor
H2_a = lambda a : H_0**2 * (Omega_r0/a**4 + (Omega_bm0+Omega_dm0)/a**3 + Omega_Lambda0 + Omega_k0/a**2)
H2_z = lambda z : H2_a(1/(1+z))
# Critical density
rho_cr_z = lambda z: 3*H2_z(z)/(8*np.pi*CGS.G)

V_r = lambda r : 4*np.pi/3 * r**3 # np.sum(R_Mpc<r)*dV
_M_r = lambda r : np.sum(dens_c[R<r])*(IA2.dL*1e3*CGS.pc)**3

def M_r(r):

    _r_list = np.arange(0,min(shape)//2)*IA2.dL*1e3*CGS.pc
    _M_r_list = np.array([_M_r(_r) for _r in _r_list])

    try:
        return RegularGridInterpolator(_r_list[np.newaxis],_M_r_list)(r)
    except:
        return RegularGridInterpolator(_r_list[np.newaxis],_M_r_list)([r])[0]
        
rhoavg_r = lambda r : M_r(r)/V_r(r)

def crop_to_centre(val):

    min_cells = tuple(min(_centre,IA2.dim-_centre) for _centre in centre)
    xlim = centre[0]-min_cells[0], centre[0]+min_cells[0]
    ylim = centre[1]-min_cells[1], centre[1]+min_cells[1]
    zlim = centre[2]-min_cells[2], centre[2]+min_cells[2]

    return val[xlim[0]:xlim[1],ylim[0]:ylim[1],zlim[0]:zlim[1]]

def get_Rmesh(shape):

    xyz = tuple(np.arange(-sh//2,sh//2,)*IA2.dL * 1e3*CGS.pc for sh in shape)
    XYZ = np.meshgrid(*xyz,indexing='ij')

    return np.sum([X**2 for X in XYZ],axis=0)**.5

if __name__=="__main__":

    from scipy.optimize import root_scalar
    from utils.enzo import IA2, IA2Data, get_redshift
    import argparse, os
    
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("Delta", help="Delta times the critical density of the universe.", type=str)
    parser.add_argument("-snapshots", help="Which snapshot to calculate parameters of. Defaults to 15th.",required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    Delta = float(args.Delta)
    snapshot = int(args.snapshot) if args.snapshot else 15

    if Path.endswith('E3A-PM'):
        IA2.dim = 1024

    z = get_redshift(snapshot)
    rho_cr = rho_cr_z(z)

    ia2 = IA2Data(Path)
    centre = ia2.get_centre(snapshot)
    rho = ia2.get_val("Density",snapshot)
    dens = rho + ia2.get_val("Dark_Matter_Density",snapshot)
    temp = ia2.get_val("Temperature",snapshot)
    prs_c = crop_to_centre(rho)/CGS.mp * CGS.kB * crop_to_centre(temp)

    dens_c = crop_to_centre(dens)
    shape = dens_c.shape
    R = get_Rmesh(shape)

    r_Delta = root_scalar(lambda r : rhoavg_r(r)/rho_cr_z(z)-Delta, x0=2*(1e6*CGS.pc),x1=2.01*(1e6*CGS.pc)).root
    M_Delta = M_r(r_Delta)
    rho_Delta = rhoavg_r(r_Delta)
    c_s = (IA2.gamma*prs_c[R<r_Delta].mean()/rho_Delta)**.5
    t_s = 2*rho_Delta/c_s

    np.savetxt(os.path.join(Path,f'cluster_{Delta}.txt'),[r_Delta,M_Delta,rho_Delta,c_s,t_s],header=f"r_{Delta} M_{Delta} rho_{Delta} c_s t_s")



