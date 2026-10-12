"""
Functions that calculate derived fields.

Functions
---------
    calc_strain, calc_stress
        Calculate the strain and stress of a fluid, respectively.

"""

import numpy as np

def calc_strain(vx,vy,vz=None,dxdydz=None,gradient_func=None):
    """
    Calculates and returns the strain given velocity fields.
    Works for 2 or 3 dimensions, assumes even spacing in every dimension.

    The strain tensor is
            \\sigma_ij = 0.5 * (\\partial_i v_j + \\partial_j v_i)

    Parameters
    ----------
        vx, vy, vz : np.ndarray
            Velocity in the x, y and z directions.
        dxdydz : iterable
            Step size between different values of field, one iterable for each dimension.
    Returns
    -------
        strain : np.ndarray
            Strain tensor.
    """

    if vz:
        ndim = 3
        v = [vx,vy,vz]
    else:
        ndim = 2
        v = [vx,vy]

    gradients = [gradient_func(_v,dxdydz) for _v in v]

    strain = []
    for i in range(ndim):
        strain.append([])
        for j in range(ndim):
            strain[i].append([])
            strain[i][j] = 0.5 * (gradients[i][j]+gradients[j][i])

    return np.array(strain)

def calc_stress(prs,visc,vx,vy,vz=None,dxdydz=None):
    """
    Calculates and returns the stress given velocity fields.
    Works for 2 or 3 dimensions, assumes even spacing in every dimension.

    Parameters
    ----------
        prs : np.ndarray
            Pressure field.
        visc : float
            Dynamic viscosity of the fluid
            (as opposed to kinematic viscosity, respectively mu and nu in Tennekes).
        vx, vy, vz : np.ndarray
            Velocity in the x, y and z directions.
        dxdydz : iterable
            Step size between different values of field, one iterable for each dimension.
    Returns
    -------
        stress : np.ndarray
            Stress tensor.
    """

    ndim = prs.ndim

    strain = calc_strain(vx,vy,vz,dxdydz)
    stress = 2*visc*strain

    for i in range(ndim):
        stress[i][i] -= prs

    return stress

def calc_ek_field(dens,vx,vy=None,vz=None):
    """ Calculate kinetic energy. """
    if not vy:
        vy = np.zeros_like(vx)
    if not vz:
        vz = np.zeros_like(vx)
    return 0.5 * dens * (vx**2 + vy**2 + vz**2)