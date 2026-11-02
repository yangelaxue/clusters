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

def extend_domain(val,shape_ext,centre):
    """ Extend domain about a centre. """

    shape_init = val.shape
    center_ext = tuple(_sh//2 for _sh in shape_ext)
    slices = []
    for d in range(len(shape_init)):
        slices.append(slice(center_ext[d]-centre[d],center_ext[d]-centre[d]+shape_init[d]))
    slices = tuple(slices)
    ret = np.zeros(shape_ext)
    ret[slices] += val
    return ret

def crop_domain(val,shape_cr,centre):
    """
    Crop a domain about a grid centre.
    If the shape_cr does not fit about the centre, just remove entire section.
    """
    if val.shape[0]<=shape_cr[0]:
        return val
    slices = []
    for cen,sh_cr in zip(centre,shape_cr):
        cr = cen-sh_cr//2
        l = max(0,cr)
        r = l+sh_cr
        slices.append(slice(l,r+1))
    slices = tuple(slices)
    val = val[slices]
    return val[slices]