"""
Calculates gravitational values relating to a density field.
Rewritten from a code from Masters project for readability and sanity

Author: Angela Xue
Date: March 2026
"""

#%% Imports.

from numpy import pi
from gradient import invlapl_FT, gradient_FT

def get_gravpot(rho, dxdydz=None, G=1.):
    """
    Calculates the gravitational potential field energy density field.
    Aside: if inputs are given in natural units, the potential is measured in units of c^2.

    Parameters
    ----------
    rho : np.ndarray
        Density field values..
    dxdydz : iterable
        Step size between different values of field, one iterable for each dimension.
    G : float, optional
        Gravitational constant

    Returns
    -------
    potential : np.ndarray
        Gravitational potential field in natural units.
    """
    
    potential = invlapl_FT(rho, dxdydz)
    potential *= 4*pi*G
    
    return potential

def get_gravforce_from_pot(potential,dxdydz):
    """
    Calculate the gravitational vector field given gravitational potential in
    cartesian coordinates by taking the gradient of said potential field.

    Parameters
    ----------
    potential : np.ndarray
        Gravitational potential field.
    dxdydz : iterable
        Step size between different values of field, one iterable for each dimension.

    Returns
    -------
    grav_force : list
        List of gravitational accceleration in the x, y and z directions respectively.
    """
    
    return [-gradient for gradient in gradient_FT(potential,dxdydz)]
    