"""
Save properties of the ENZO Itasca IA2 clusters.

Parameters saved are:
    Cluster center of mass index.
    Cluster parameters (r_Delta, M_Delta, c_s, t_s).
    Cluster centre density.

Author: Angela Xue
Date:
Revised: October 2026
"""

import numpy as np
import os
from utils.units import CGS

def save_centerofmass_indices(ia2):
    """ Calculate centre of mass indices. """

    # from utils.clusters import calc_centreofmass
    from scipy.ndimage import center_of_mass
    rho_gen = ia2.gen_vals("Density",snapshots=snapshots,c=c)
    dm_gen = ia2.gen_vals("Dark_Matter_Density",snapshots=snapshots,units="code",c=c)
    line = "{:02d} {:4.1f} {:4.1f} {:4.1f}\n"
    with open(os.path.join(ia2.Path,"centreindices.txt"),"w") as f:
        f.write("s  idx   idy   idz\n")
        for s, rho, dm in zip(snapshots,rho_gen,dm_gen):
            dens = rho + dm
            centre = center_of_mass(dens)
            # centre = calc_centreofmass(dens)
            # print(line.format(s,*centre))
            f.write(line.format(s,*centre))

def save_centredensity(ia2):
    """ Save cluster centre density values. """

    line = "{:02d} {:8.7e}\n"
    with open(os.path.join(ia2.Path,"centredensity.txt"),"w") as f:
        f.write("s  density\n")
        for s in snapshots:
            centre = ia2.get_centre(snapshot=s,c=c)
            centre = tuple(int(cen) for cen in centre)
            density = ia2.get_val("Density",snapshot=s,units="cgs",c=c)[centre]
            density += ia2.get_val("Dark_Matter_Density",snapshot=s,units="cgs",c=c)[centre]
            f.write(line.format(s,density))

def save_clusterunits(ia2):
    """ Cluster parameters (r_Delta, M_Delta, c_s, t_s) in CGS units. """

    from utils.grid import get_Rmesh
    from scipy.interpolate import RegularGridInterpolator
    from scipy.optimize import root_scalar

    Deltas = [100,200,500,1000]

    # Define critical density as a function of redshift.
    H2_a = lambda a : IA2.H_0**2 * (IA2.Omega_r0/a**4 + (IA2.Omega_bm0+IA2.Omega_dm0)/a**3 + IA2.Omega_Lambda0 + IA2.Omega_k0/a**2)
    H2_z = lambda z : H2_a(1/(1+z))
    rho_cr_z = lambda z: 3*H2_z(z)/(8*np.pi*CGS.G)
    z = IA2.get_redshift(s=15)
    rho_cr = rho_cr_z(z)

    # Define rhoavg_r
    Mpc = 1e6*CGS.pc
    dL = IA2.dL/(1+z)*c / Mpc # Physical units, Mpc
    dens = ia2.get_val("Density",snapshot=15,c=c,units='cgs')
    dens += ia2.get_val("Dark_Matter_Density",snapshot=15,c=c,units='cgs')
    dens /= rho_cr # Units of critical density.
    R = get_Rmesh((ia2.dim//c,)*3,ia2.get_centre(snapshot=15,c=c))*dL
    r_arr = np.arange(dL,R.max(),dL)
    M_arr = np.array([np.sum(dens[R<_r]) for _r in r_arr]) * dL**3 # Mass in a Mpc volume
    def M_r(r):
        try:
            return RegularGridInterpolator(r_arr[np.newaxis],M_arr)(r)
        except:
            return RegularGridInterpolator(r_arr[np.newaxis],M_arr)([r])[0]
    V_r = lambda r : 4*np.pi/3 * r**3
    rhoavg_r = lambda r : M_r(r) / V_r(r)
    # Calculate cluster units in CGS.
    line = "{:5d} {:8.7e} {:8.7e} {:8.7e} {:8.7e}\n"
    with open(os.path.join(ia2.Path,"clusterunits.txt"),"w") as f:
        f.write("Delta r_Delta       M_Delta       c_s           t_s\n")
        for Delta in Deltas:
            r_Delta = root_scalar(lambda r : rhoavg_r(r)-Delta,x0=1).root # Mpc
            M_Delta = M_r(r_Delta) # rho_cr in 1 Mpc volume
            M_Delta *= rho_cr * Mpc**3 # g
            temp = np.array(ia2.get_val("Temperature",snapshot=15,c=c,units='cgs'),dtype=np.float64)
            c_s = np.mean((5/3*CGS.kB*temp[R<r_Delta]/CGS.mp)**.5) # cm/s
            r_Delta *= Mpc # cm
            t_s = 2*r_Delta/c_s # s
            f.write(line.format(Delta,r_Delta,M_Delta,c_s,t_s))

if __name__=="__main__":

    from utils.ia2 import IA2, IA2Data

    # Parse arguments.
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("cluster", help="Which ENZO Itasca IA2 cluster to save the properties for.",type=str)
    parser.add_argument("-DEBUG", help="Switch to debug mode.",type=bool)
    parser.add_argument("-c",type=int)
    args = parser.parse_args()
    cluster = args.cluster
    DEBUG = args.DEBUG

    ia2 = IA2Data(cluster)
    c = 1
    snapshots = IA2.snapshots[ia2.cluster]

    if DEBUG:
        c = args.c if args.c else 10
        snapshots = IA2.snapshots[ia2.cluster]

    save_centerofmass_indices(ia2)
    # save_centredensity(ia2)
    # save_clusterunits(ia2)