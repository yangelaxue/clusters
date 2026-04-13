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

def calc_avgvelocity(dens,vxvyvz):

    vel = np.sum(dens*np.array(vxvyvz),axis=(1,2,3))/np.sum(dens)

    return vel