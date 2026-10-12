"""
Functions that deal with generation or manipulations of grids.

Author: Angela Xue
Date: October 2026
"""

import numpy as np

def get_xyz(shape,centre=None):
    if not centre:
        centre = tuple(sh//2 for sh in shape)
    xyz = tuple(np.arange(0,sh)-cen for sh,cen in zip(shape,centre))
    return xyz

def get_meshgrid(shape,centre=None):
    if not centre:
        centre = tuple(sh//2 for sh in shape)
    xyz = tuple(np.arange(0,sh)-cen for sh,cen in zip(shape,centre))
    XYZ = np.meshgrid(*xyz,indexing='ij')

def get_Rmesh(shape,centre=None):
    """ Generate a spherical radial meshgrid. """
    if not centre:
        centre = tuple(sh//2 for sh in shape)
    xyz = tuple(np.arange(0,sh)-cen for sh,cen in zip(shape,centre))
    XYZ = np.meshgrid(*xyz,indexing='ij')
    return np.sum([X**2 for X in XYZ],axis=0)**.5