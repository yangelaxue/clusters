"""
Integrate raw data along a line of sight to calculate either the X-ray surface brightness or Compton-y SZ surface brightness.

Functions
---------

    get_XraySB()

    get_SZSB()

Author: Angela Xue
Date: Oct 2025
Revised: Oct 2026 -- combined dependent functions.
"""

import numpy as np
import os
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson

from utils.units import CGS

def get_CoolingFunction(Z):
    """
    Return cooling function.

    Metallicity, in units of solar metallicity.
    Z = 0.3, 1.0, 3.0
    """

    cooltab_fname = os.path.join(os.environ['PLUTO_DIR'],'Src/Cooling/TABULATED/cooltab_z03z3.dat')
    if float(Z)==0.3:
        T_arr, cool_arr = np.loadtxt(cooltab_fname,usecols=(0,1),skiprows=1)[:-1].T
    elif float(Z)==1.0:
        T_arr, cool_arr = np.loadtxt(cooltab_fname,usecols=(0,2),skiprows=1)[:-1].T
    elif float(Z)==3.0:
        T_arr, cool_arr = np.loadtxt(cooltab_fname,usecols=(0,3),skiprows=1)[:-1].T
    else:
        ValueError("Z must be 0.3, 1.0 or 3.0.")
    
    return CubicSpline(T_arr, cool_arr)

def get_XraySB(rho,prs=None,T=None,mu=1.,Z=0.3,dx=1.,los='x'):
    """
    Given field values and line-of-sight, integrate to calculate the X-ray surface brightness.
    Gives the CGS surface brightness value if inputs are in CGS units.

    Parameters
    ----------
    rho : np.ndarray
        Density field.
    prs, T : np.ndarray
        Pressure and temperature fields, respectively.
        Must provide at least one of these.
    mu : float
        Mean molecular mass.
    Z : float.
        Metallicity, in units of solar metallicity.
    dx : float
        Unit of integration.
    los : int
        Line-of-sight to integrate, los=0,1,2.
    """

    if not T:
        try:
            T = prs/(rho/CGS.mp/mu*CGS.kB)
        except:
            raise ValueError("Must provide either pressure or temperature values.")

    CoolingFunction = get_CoolingFunction(Z)

    if los in {'x','X'}:
        los = 0
    elif los in {'y','Y'}:
        los = 1
    elif los in {'z','Z'}:
        los = 2
    else:
        ValueError("Specify the line of sight as 0, x, X, 1, y, Y, 2, z or")

    SB = simpson((rho/CGS.mp/mu)**2*CoolingFunction(T),axis=los,dx=dx)

    return SB

def get_SZSB(prs,dx=1.,los='x'):
    """
    Given field values and line-of-sight, integrate to calculate the X-ray surface brightness.
    Gives the CGS surface brightness value if inputs are in CGS units.

    Parameters
    ----------
    prs : np.ndarray
        Pressure fields.
        Must be given in physical units.
    dx : float
        Unit of integration.
    los : int
        Line-of-sight to integrate, los=0,1,2.
    """

    SB = CGS.sigmaT/(CGS.me*CGS.c**2) * simpson(prs,axis=los,dx=dx)

    return SB