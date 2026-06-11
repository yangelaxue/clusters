"""
Functions that calculate important cluster properties regardless of data format.

Author: Angela Xue
Date April 2026
"""

import numpy as np

def calc_centreofmass(dens,XYZ=None,):

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

    vxvyvz = [vx]
    axis = [1]
    if vy:
        vxvyvz.append(vy)
        axis.append(2)
    if vz:
        vxvyvz.append(vz)
        axis.append(3)
    vxvyvz = tuple(vxvyvz)
    axis = tuple(axis) if len(axis)>1 else None

    vel = np.array([np.sum(dens*vx,axis=axis) for vx in vxvyvz])/np.sum(dens)

    return vel