"""
Functions that calculate important cluster properties regardless of data format.

Functions
---------
    calc_centreofmass
        Calculate the centre of mass position of a cluster.
    calc_avgvelocity
        Calculate the mass-weighted average velocity of a cluster.

Author: Angela Xue
Date: April 2026
Revised: October 2026
"""

import numpy as np

def get_M_r(dens,centre,dxdydz:float=None):
    """ Calculate the mass enclosed within a radius r. """
    from scipy.interpolate import RegularGridInterpolator
    from utils.grid import get_Rmesh
    if not dxdydz:
        dxdydz = 1.
    R = get_Rmesh(dens.shape,centre)*dxdydz
    dr = (dxdydz**2*3)**.5
    dV = dxdydz**3
    r_arr = np.arange(dr,R.max(),dr)
    M_arr = np.array([np.sum(dens[R<_r]) for _r in r_arr]) * dV # Mass in a Mpc volume
    return RegularGridInterpolator(r_arr[np.newaxis],M_arr)

def NFW_profile(r,n_0,r_s):
    """ Navarro-Frenk-White profile. """
    ret = np.zeros_like(r)
    zeros = np.where(r==0)
    non_zeros = np.where(r!=0)
    ret[zeros] = n_0
    ret[non_zeros] = n_0 / (r[non_zeros]/r_s) / (1+r[non_zeros]/r_s)**2
    return ret

def gNFW_profile(r,n_0,r_500,c_500,alpha,beta,gamma):
    """ generalised Navarro-Frenk-White profile. """
    r_s = r_500/c_500
    x = r/r_s
    ret = np.zeros_like(r)
    zeros = np.where(r==0)
    non_zeros = np.where(r!=0)
    ret[zeros] = n_0
    ret[non_zeros] = n_0 / (c_500*x[non_zeros])**(gamma) / (1+(c_500*x[non_zeros])**alpha)**((beta-gamma)/alpha)
    return ret

def beta_profile(r,n_0,r_c,beta):
    """ beta-profile of a cluster. """
    return n_0 * (1+(r/r_c)**2)**(-3*beta/2)

def einasto_profile(r,n_2,r_2,alpha):
    """ Einasto profile. """
    return n_2 * np.exp(-2/alpha * ((r/r_2)**alpha-1))
    
def calc_centreofmass(dens):
    """ Calculate the centre of mass position of a cluster. """
    from scipy.ndimage import center_of_mass
    return center_of_mass(dens)

def calc_avgvelocity(dens,vx,vy=None,vz=None):
    """
    Calculate the mass-weighted average velocity.

    Parameters
    ----------
    dens : np.ndarray
        Density or mass field.
    vxvyvz : tuple
        Tuple of velocities. Can be however many dimensions.
    """
    v = [vx]
    axis = [1]
    if vy:
        v.append(vy)
        axis.append(2)
    if vz:
        v.append(vz)
        axis.append(3)
    v = tuple(v)
    axis = tuple(axis) if len(axis)>1 else None
    vel = np.array([np.sum(dens*_v,axis=axis) for _v in v])/np.sum(dens)
    return vel