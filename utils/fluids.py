"""
Calculate and return the rms values.
Intended to be agnostic toward simulation programme.

Functions
---------
    get_mach_rms
        Takes all raw data given by any simulation to calculate turbulent Mach number for a single timeslice.

    get_c_s
        Calculate sound speed in an adiabatic fluid.

    get_val_rms
        Calculate rms value for a given variable for a single timeslice.
    

Author: Angela Xue
Date: March 2026
"""

import numpy as np
from scipy.integrate import simpson
from gradient import gradient_discrete, gradient_FT

def get_mach_rms(vx1,vx2=None,vx3=None,subtract_mean=False,normalise=0,**kwargs):
    """
    Takes all raw data given by any simulation to calculate turbulent Mach number for a single timeslice.

    Parameters
    ----------
    vx1, vx2, vx3 : np.ndarray
        Velocity in the x, y and z directions.
    subtract_mean : bool
        True to use delta v = v-v_mean in calculations.
        False to use v in calculations.
    normalise : int
        Whether to calculate the root-mean-sqaure Mach number or not.
        0 to not calculate the Mach number.
        1 to calculate mach_rms = v_rms/c_s_mean.
        2 to calculate mach_rms = (v/c_s)_rms.
    
    Returns
    -------
    v_rms : float
    OR
    mach_rms & c_s : float
    """

    # Ensure all arrays are of the same shape.
    if type(vx2)==np.ndarray: assert vx1.shape==vx2.shape, "Shapes of vx1 and vx2 must match."
    else:
        vx2 = np.zeros_like(vx1)
        vx3 = np.zeros_like(vx1)
    if type(vx3)==np.ndarray: assert vx1.shape==vx3.shape, "Shapes of vx1 and vx3 must match."
    else: 
        vx3 = np.zeros_like(vx1)

    # Calculate the velocity magnitude field
    if subtract_mean:
        v = ((vx1-vx1.mean())**2 + (vx2-vx2.mean())**2 + (vx3-vx3.mean())**2)**.5
    else:
        v = (vx1**2 + vx2**2 + vx3**2)**.5

    # normalise
    if normalise==0:
        v_rms = np.mean(v**2)**.5
        return v_rms
    else:
        # Ensure that rho and prs are defined.
        try:
            rho = kwargs['rho']
            prs = kwargs['prs']
            gamma = kwargs['gamma'] if 'gamma' in kwargs.keys() else 5/3 # Default to 5/3
        except:
            raise KeyError("Must provide rho and prs value to normalise against sound speed.")
        
        assert vx1.shape==rho.shape, "Shapes of rho and vx1 must match."
        assert vx1.shape==prs.shape, "Shapes of prs and vx1 must match."
        
        if normalise==1:
            c_s = get_c_s(rho,prs,gamma)
            mach_rms = np.mean(v**2)**.5/c_s.mean()
        elif normalise==2:
            c_s = get_c_s(rho,prs,gamma)
            mach_rms = np.mean(v**2/c_s**2)**.5
        else:
            raise NotImplementedError(f"Parameter normalise must be set to 0, 1 or 2.")
        return mach_rms, c_s.mean()

def get_c_s(rho,prs,gamma=5/3):
    """ Calculate average adiabatic sound speed in a field. """
    c_s = (gamma*prs/rho)**.5
    return c_s

def get_val_rms(val,subtract_mean=False,normalise=0):
    """
    Calculate rms value for a given variable for a single timeslice.

    Parameters
    ----------
    val : np.ndarray
        Field values to take the root-mean-square of.
    subtract_mean : bool
        True to use delta val = val-val_mean in calculations.
        False to use val in calculations.
    normalise : int
        0 to just calculate the root-mean-square.
        1 to divide the root-mean-square value by val.mean().
        2 to divide the field values by its mean before taking the root-mean-square.
    
    Returns
    -------
    val_rms : float
    """

    val_mean = val.mean()

    # subtract_mean
    if subtract_mean:
        val -= val_mean

    # normalise
    if normalise==0:
        val_rms = np.mean(val**2)**.5
    elif normalise==1:
        val_rms = np.mean(val**2)**.5/val_mean
    else:
        NotImplementedError(f"Parameter normalise must be set to 0 or 1.")

    return val_rms

def get_strain(vx1,vx2,vx3,dxdydz=None,gradient_func=gradient_FT):
    """
    Calculates and returns the strain given velocity fields.
    Works for 2 or 3 dimensions, assumes even spacing in every dimension.

    The strain tensor is
            \\sigma_ij = 0.5 * (\\partial_i v_j + \\partial_j v_i)

    Parameters
    ----------
        vx1, vx2, vx3 : np.ndarray
            Velocity in the x, y and z directions.
        dxdydz : iterable
            Step size between different values of field, one iterable for each dimension.
        
    Returns
    -------
        strain : np.ndarray
            Strain tensor.
    """

    if vx3:
        ndim = 3
        assert vx1.shape==vx2.shape and vx1.shape==vx3.shape, "Velocity fields must be of the same shape"
        v = [vx1,vx2,vx3]
    else:
        ndism = 2
        assert vx1.shape==vx2.shape, "Velocity fields must be of the same shape"
        v = [vx1,vx2]

    gradients = [gradient_func(vx,dxdydz) for vx in v]

    strain = []
    for i in range(ndim):
        strain.append([])
        for j in range(ndim):
            strain[i].append([])
            strain[i][j] = 0.5 * (gradients[i][j]+gradients[j][i])

    return np.array(strain)


def get_stress(prs,visc,vx1,vx2,vx3=None,dxdydz=None):
    """
    Calculates and returns the stress given velocity fields.
    Works for 2 or 3 dimensions, assumes even spacing in every dimension.

    Parameters
    ----------
        prs : np.ndarray
            Pressure field.
        visc : float
            Dynamic viscosity of the fluid (as opposed to kinematic viscosity,
                                            respectively mu and nu in Tennekes).
        vx1, vx2, vx3 : np.ndarray
            Velocity in the x, y and z directions.
        dxdydz : iterable
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
    stress = 2*visc*strain.copy()

    for i in range(ndim):
        stress[i][i] -= prs

    return stress

def get_timeaverage(vals,times=None):
    """
    Calculates the time-averaged values for a timeseries.
    If times are not provided, even time spacing between snapshots will be assumed.

    Parameters
    ----------
        vals : iterable
            Time series of field values to average.
            Must have shape (ntimes, *val.shape)
        times : iterable
            Times that correlate to the snapshots of vals
    
    Return
    ------
        np.ndarray
            A time-averaged field.
    """

    if not times:
        times = np.arange(len(vals))
    T = times[-1]-times[0]
    return 1/T * simpson(vals,times,axis=0)

def get_correlation(val1,val2,times=None):
    """ Calculates the correlation between two variables. """

    assert val1.shape==val2.shape, "Fields must be of the same shape."

    num = get_timeaverage(val1*val2,times)

    val1_avg = get_timeaverage(val1**2,times)
    val2_avg = get_timeaverage(val2**2,times)
    den = (val1_avg*val2_avg)**.5

    return num/den

def get_probabilitydensityfunction(val,weight,inc=0,axis=0,normalise=False):
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
            number of 
        inc : int
            If pdf of increments, then subtract the arrays by some increment.
        axis : int
            If pdf of increments, determine the asix which to take increments.
        normalise : bool
            Normalise values in val to its room-mean-square.
    """

    if inc:
        val -= np.roll(val,inc,axis=axis)
    if normalise:
        val /= (np.mean((val-val.mean())**2))**.5

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

"""
    # if bins:
    #     # If bins are provided.
    #     if type(bins)==int:
    #         bins = np.linspace(val.min(),val.max(),bins,endpoint=False)
    #     else:
    #         assert len(bins)>1, "Provided bins must be arraylike."
    #     nbins = len(bins)
    #     pdf = np.zeros(nbins)

        for i in range(nbins-1):
            pdf[i] = np.logical_and(bins[i]<=val,val<bins[i+1]).sum() / (bins[i+1]-bins[i])
        pdf[-1] = ((bins[-1]<=val)).sum() / (bins[-1]-bins[-2])
        pdf /= val.size

    else:
        # Else, weights are used.
        nbins = val.size // weight
        bins = np.zeros(nbins)

        val_sorted = np.sort(val.flatten())
        pdf = np.zeros(nbins)
        for i in range(nbins):
            idx = i * weight
            bins[i] = np.mean(val_sorted[idx:idx+weight])
            pdf[i] = weight / (val_sorted[idx+weight-1]-val_sorted[idx])
        pdf /= val.size

    return bins, pdf
"""



from functools import reduce
import operator

def calc_pdf(ar,min=99999,max=99999,weight=100,inc=0,ax=0,Normalized=False):
    if len(ar) == 0:
      print('No array provided! Exiting!')
      return
    if min == 99999:
      min=ar.min()
    if max == 99999:
      max=ar.max()
   # If PDF of increment, then increment the array
    if inc:
      ar -= np.roll(ar,inc,axis=ax)
#    # Find the total length of data set
    arsize=ar.size
#    # Find the RMS value of data set and normalize to it.
    if Normalized == True:
      rmsval = np.sqrt(np.mean((ar-np.mean(ar))**2))
      if rmsval != 0:
         ar = ar/rmsval
#    # Reshape the array to 1D & sort it.
    arr=ar.flatten()
    print(arr)
    arr = np.sort(arr,kind='heapsort')
    print(arr)
#    # Empty arrays for output
#    bins=int(arsize/weight); pdf=np.zeros(bins); binvals=np.zeros(bins)
#    # Fill the bins 
#    for i in range(bins):
#       start=i*weight
#       binvals[i] = np.mean(arr[start:start+weight])
#       pdf[i] = weight/(arr[start:start+weight].max()-arr[start:start+weight].min())
#    pdf = pdf/arsize
#    return binvals,pdf

def main():

    # val = np.random.random((128,128,128))
    val = np.array([8,3,5,1,7,2,4,0,6]).reshape((3,3))

    # calc_pdf(val)
    get_probabilitydensityfunction(val,bins=[])

if __name__=='__main__':
    main()

