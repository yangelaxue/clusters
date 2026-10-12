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

def calc_centreofmass(dens,XYZ=None,):
    """ Calculate the centre of mass position of a cluster. """
    if not XYZ:
        xyz = tuple(np.arange(0,sh) for sh in dens.shape)
        XYZ = np.meshgrid(*xyz,indexing='ij')
    centre = np.sum(dens*np.array(XYZ),axis=(1,2,3))/np.sum(dens)
    return centre

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