"""
Key turbulence analysis.

Functions
---------
    calc_c_s
        Calculate the adiabatic sound speed of a fluid.
        Spatially dependent.
    calc_mach_rms
        Calculate turbulent Mach number.
        Can return v_rms, v_rms/c_s_rms or (v/c_s)_rms depending on c_s input.
    calc_val_rms
        Calculate rms of any field values.
    calc_probabilitydensityfunction
        Calculate the probability density function of field values.

For statistically static turbulence...
    calc_timeaverage
        Calculate the time-average field values.
    calc_correlation
        Calculate the time correlation between two field values.

Author: Angela Xue
Date: March 2026
Revised: October 2026
"""

import numpy as np

def calc_c_s(rho,prs,gamma=5./3.):
    """ Calculate average adiabatic sound speed in a field. """
    c_s = (gamma*prs/rho)**.5
    return c_s

def calc_mach_rms(vx,vy=None,vz=None,c_s=1.):
    """
    Takes all raw data given by any simulation to calculate turbulent Mach number for a single timeslice.

    Parameters
    ----------
        vx, vy, vz : np.ndarray
            Velocity in the x, y and z directions.
        c_s : float, np.ndarray
            Sound speed.
            If c_s==1., v_rms is returned.
            If type(c_s)==float, v_rms/c_s_rms is returned.
            If type(c_s)==np.ndarray, (v/c_s)_rms is returned.
    Returns
    -------
    v_rms : float
    """

    if not vy:
        vy = np.zeros_like(vx)
    if not vz:
        vz = np.zeros_like(vx)

    v = (vx**2 + vy**2 + vz**2)**.5
    v_rms = np.mean(v**2/c_s**2)**.5

    return v_rms

def calc_val_rms(val):
    """ Calculate rms value for a given variable for a single timeslice. """
    return np.mean(val**2)**.5

def calc_probabilitydensityfunction(val,weight,increment,axis):
    """
    Calculate the probability density function of a field.
        The percentage of points expected to fall within x and x+dx is P(x)dx
        If we weight this by 100 data points per bin, then
                weight = P(x)dx
        from which we can calculate
                P(x) = weight/dx
        where dx is the distance between data points within this bin.
    Based on Tulasi's function calc_pdf from
    https://github.com/tulasinandan/TurbAn/blob/master/Analysis/Simulations/AnalysisFunctions.py

    Parameters
    ----------
        val : np.ndarray
        weight : int
            Number of elements to be in a single bin in the PDF.
        inc : int
            If PDF of increments, then subtract the arrays by some increment.
        axis : int
            If PDF of increments, determine the asix which to take increments.
    """

    if increment:
        val -= np.roll(val,increment,axis=axis)

    val_sorted = np.sort(val.flatten())
    nbins = val.size//weight
    pdf = np.zeros(nbins)
    binvals = np.zeros(nbins)

    for i in range(nbins):
        start_idx = i*weight
        binvals[i] = np.mean(val_sorted[start_idx:start_idx+weight])
        pdf[i] = weight / (val_sorted[start_idx:start_idx+weight].max() - val_sorted[start_idx:start_idx+weight].min())
    pdf /= val.size

    return binvals, pdf

#%% Below for statistically stable turbulence.

def calc_timeaverage(vals,times=None):
    """
    Calculates the time-averaged values for a timeseries.

    Parameters
    ----------
        vals : iterable
            Time series of field values to average.
            Must have shape (ntimes, *val.shape)
        times : iterable
            Times that correlate to the snapshots of vals.
            If not provided, times will be spaced apart by 1.
    
    Return
    ------
        np.ndarray
            A time-averaged field.
    """
    from scipy.integrate import simpson

    if not times:
        times = np.arange(len(vals))
    T = times[-1]-times[0]
    return 1/T * simpson(vals,times,axis=0)

def calc_correlation(val1,val2,times=None):
    """ Calculates the correlation between two variables. """

    assert val1.shape==val2.shape, "Fields must be of the same shape."

    numerator = calc_timeaverage(val1*val2,times)

    val1_avg = calc_timeaverage(val1**2,times)
    val2_avg = calc_timeaverage(val2**2,times)
    denomenator = (val1_avg*val2_avg)**.5

    return numerator/denomenator