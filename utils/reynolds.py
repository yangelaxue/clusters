"""
Calculation of basic turbulence numbers.
Based on A First Course in Turbulence, by H. Tennekes and J. L. Lumley, chapter 2.

Reynold's Decomposition
-----------------------
Fluctuating values are decomposed into their time-averaged mean and their fluctuation:
    \\tilda u = \\tilda u(t,x) = U(x) + u(t,x) = U + u
    U = lim_{T -> \\infinity} 1/T * \\int_0^T \\tilda u dt

------------------
Author: Angela Xue
"""

#%% Imports.

import numpy as np
from scipy import integrate

from utils.differentiation_utils import gradient_FT

#%% Reynold's decomposition.

def get_timeaverage(vals, times):
    """
    Calculates the time-averaged values for a timeseries.
    
    Parameters
    ----------
        vals : array_like, shape==(len(times), <grid_shape>)
            A time ordered arrays of all the grid values of a variable.
        times : array_like
            An array of all the time values corresponding to vals.
    Returns
    -------
        array_like
            A time-averaged field with shape==vals[0].shape.
    """

    T = times[-1]-times[0]
    return 1/T * integrate.simpson(vals,x=times,axis=0)

def get_fluctuations(vals,times,t0=None):
    """
    Calculates the fluctuating values for a timeseries.
    
    Parameters
    ----------
        vals : array_like, shape==(len(times), <grid_shape>)
            A time ordered arrays of all the grid values of a variable.
        times : array_like
            An array of all the time values corresponding to vals.
        t0 : float
            A time in the range of times which to begin the calculation of the time-average velovity.
    Returns
    -------
        array_like
            An afieldrray of fluctuations with shape==vals.shape.
    """
    t0 = times[0] if not t0 else t0

    timeaverage = get_timeaverage(vals[times>=t0],times[times>=t0])
    vals_fluctuations = np.array([val-timeaverage for val in vals])

    return vals_fluctuations

def get_strain(vx1,vx2,vx3=None,dxdydz=None):
    """
    Calculates and returns the strain given velocity fields.
    Works for 2 or 3 dimensions, assumes even spacing in every dimension.

    Parameters
    ----------
        vx1, vx2, vx3 : np.ndarray
            Velocity fields, vx3 is optional.
        dxdydz : tuple
            Step size between different values of field, one iterable for each dimension.

    Returns
    -------
        strain : np.ndarray
            Strain tensor.
    """

    if vx3:
        ndim = 3
        assert vx1.shape==vx2.shape and vx1.shape==vx3.shape, "velocity fields must be of the same shape"
        v = [vx1,vx2,vx3]
    else:
        ndim = 2
        assert vx1.shape==vx2.shape, "velocity fields must be of the same shape"
        v = [vx1,vx2]
    if not dxdydz:
        dxdydz = tuple(1. for _ in range(ndim))
    
    strain = [[0 for _i in range(ndim)] for _j in range(ndim)]
    gradients = []

    for _v in v:
        grad = gradient_FT(_v,dxdydz)
        gradients.append(grad)
        
    for i in range(ndim):
        for j in range(ndim):
            strain[i][j] = 0.5 * (gradients[i][j]+gradients[j][i])
    
    return np.array(strain)

def get_stress(prs,viscosity,vx1,vx2,vx3=None,dxdydz=None):
    """
    Calculates and returns the stress using 'get_strain'.
    Works for 2 or 3 dimensions, assumes even spacing in every dimension.

    Parameters
    ----------
        prs : np.ndarray
            Pressure field.
        viscosity : float
            Dynamic viscosity of the fluid (as opposed to kinematic viscosity,
                                            respectively mu and nu in Tennekes).
        vx1, vx2, vx3 : np.ndarray
            Velocity fields, vx3 is optional.
        dxdydz : tuple
            Step size between different values of field, one iterable for each dimension.

    Returns
    -------
        stress : np.ndarray
            Stress tensor.
    """

    if prs.ndim==3:
        assert vx3, "3 dimensions require vx3"
        ndim = 3
    else:
        ndim = 2

    strain = get_strain(vx1,vx2,vx3,dxdydz)
    stress = 2*viscosity*strain.copy()

    for i in range(ndim):
        stress[i][i] -= prs

    return stress

def get_correlation(val1,val2,times):
    """
    Calculates the correlation between two variables.
    
    Parameters
    ----------
        val1, val2 : array_like, shape==(len(times), <grid_shape>)
            Time ordered arrays of all the grid values of a variable.
        times : array_like
            An array of all the time values corresponding to vals.
    Returns
    -------
        array_like
            A correlation field with shape==vals[0].shape.
    """

    assert val1.shape==val2.shape, "Fields must be of the same shape."

    num = get_timeaverage(val1*val2,times)

    val1_avg = get_timeaverage(val1**2,times)
    val2_avg = get_timeaverage(val2**2,times)
    den = (val1_avg*val2_avg)**.5

    return num/den


