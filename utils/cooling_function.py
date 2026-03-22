"""
Define the cooling rate function from "$PLUTO_DIR/Src/Cooing/TABULATED/cooltable.dat".

Date: March 2026
"""

import numpy as np
import os
from scipy.interpolate import CubicSpline

def get_CoolingFunction(Z=0.3):
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