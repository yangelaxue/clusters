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