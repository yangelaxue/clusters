"""
For calculating the probability densities and related functions.
Based on A First Course in Turbulence, by H. Tennekes and J. L. Lumley, chapter 6.

The probability density function is given by:
    B(u)\\du = lim_{T -> \\infinity} 1/T * \\Sum \\Delta t

Author: Angela Xue
"""

#%% Imports

import numpy as np

#%% Begin

def get_pdf(vals,times=None,nbins=10):
    """
    Get the probability density function of a field of scalars.

    Parameters
    ----------
        vals : np.ndarray
            A time series of a scalar field.
        times : arraylike
            A list or array of times correspoinding to vals.
        nbins : int
            Number of bins which to evaluate the probability function.
    
    Returns
    -------
    """

    if not times:
        times = np.arange(len(vals))

    bins = np.linspace(vals.min(),vals.max(),nbins)

    pdf = np.zeros((nbins,*vals.shape))

    for i, bin in enumerate(bins):
        for j,val in enumerate(vals):
            where = np.where(bin<=val and val<bins[i+1])
            pdf[i][where] += times[i]
    
    pdf /= times[-1]-times[0]

    return bins, pdf
                


