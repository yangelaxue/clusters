"""
Functions to calculate gradients either using stencils or Fourier transforms.
Rewritten from a November 2022 code from Masters project for readability.

Functions
---------
    gradient_discrete
        Calculates the gradient of a field with stencils.
    lapl_discrete
        Applies the Laplace operator to a field with a standard stencil.
    gradient_FT
        Calculates the gradient of a field with Fourier and inverse transform.
    invlapl_FT
        Calculates the inverse Laplace transform os a field with Fourier and inverse transform.

Author: Angela Xue
Date: March 2026
"""

import numpy as np
from scipy.fftpack import fftn, ifftn

def gradient_discrete(field, dxdydz=None, stencil:int=3):
    """
    Calculates the gradient of a field of any dimensions using a 3, 5, 7 or 9 point
        stencil. This function assumes periodic boundary conditions and even spacing
        between points.

    Parameters
    ----------
    field : np.ndarray
        Field values in arbitrary dimentions.
    dxdydz : iterable
        Step size between different values of field, one iterable for each dimension.
    stencil : int, optional
        Type of stencil used to compute the gradient. The default is 3.

    Returns
    -------
    gradient : list
        List of gradients in the different directions of field.
    """

    # If dxdydz is not provided, set to tuple of ones.
    try:
        assert len(dxdydz)==field.ndim
    except:
        dxdydz = (1,)*field.ndim
    
    f = field.copy()
    
    # Assume periodic boundary conditions by pasting
    # the edge faces on opposite sides of the volume.
    for axis, dimension in enumerate(f.shape):
        
        f_slice, l_slice = (
            [slice(0,dim) for dim in f.shape],
            [slice(0,dim) for dim in f.shape],)
        
        f_slice[axis] = slice(0,stencil//2)
        l_slice[axis] = slice(dimension-stencil//2,dimension)
    
        f = np.concatenate((f[tuple(l_slice)], f, f[tuple(f_slice)]), axis=axis)
        
    gradient = []
    for axis, (dx,dimension) in enumerate(zip(dxdydz,f.shape)):
        
        if stencil==3:
            f_1, l_1 = (
                [slice(1,-1) for d in f.shape],
                [slice(1,-1) for d in f.shape],
            )
            f_1[axis] = slice(0,-2)
            l_1[axis] = slice(2,dimension)
            gradient.append((f[tuple(l_1)]-f[tuple(f_1)])/(2*dx))
            
        elif stencil==5:
            f_1, f_2, l_1, l_2 = (
                [slice(2,-2) for d in f.shape],
                [slice(2,-2) for d in f.shape],
                [slice(2,-2) for d in f.shape],
                [slice(2,-2) for d in f.shape],
            )
            f_1[axis] = slice(1,-3)
            f_2[axis] = slice(0,-4)
            l_1[axis] = slice(3,dimension-1)
            l_2[axis] = slice(4,dimension)
            gradient.append(
                (
                    (f[tuple(l_1)]-f[tuple(f_1)])*2/3
                    + (-f[tuple(l_2)]+f[tuple(f_2)])*1/12
                )/dx
            )
            
        elif stencil==7:
            
            f_1, f_2, f_3, l_1, l_2, l_3 = (
                [slice(3,-3) for d in f.shape],
                [slice(3,-3) for d in f.shape],
                [slice(3,-3) for d in f.shape],
                [slice(3,-3) for d in f.shape],
                [slice(3,-3) for d in f.shape],
                [slice(3,-3) for d in f.shape],
            )
            f_1[axis] = slice(2,-4)
            f_2[axis] = slice(1,-5)
            f_3[axis] = slice(0,-6)
            l_1[axis] = slice(4,dimension-2)
            l_2[axis] = slice(5,dimension-1)
            l_3[axis] = slice(6,dimension)
            gradient.append(
                (
                    (f[tuple(l_1)]-f[tuple(f_1)])*3/4
                    + (-f[tuple(l_2)]+f[tuple(f_2)])*3/20
                    + (f[tuple(l_3)]-f[tuple(f_3)])*1/60
                )/dx
            )
            
        elif stencil==9:
            f_1, f_2, f_3, f_4, l_1, l_2, l_3, l_4 = (
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
                [slice(4,-4) for d in f.shape],
            )
            f_1[axis] = slice(3,-5)
            f_2[axis] = slice(2,-6)
            f_3[axis] = slice(1,-7)
            f_4[axis] = slice(0,-8)
            l_1[axis] = slice(5,dimension-3)
            l_2[axis] = slice(6,dimension-2)
            l_3[axis] = slice(7,dimension-1)
            l_4[axis] = slice(8,dimension)
            gradient.append(
                (
                    (f[tuple(l_1)]-f[tuple(f_1)])*4/5
                    + (-f[tuple(l_2)]+f[tuple(f_2)])*1/5
                    + (f[tuple(l_3)]-f[tuple(f_3)])*4/105
                    + (-f[tuple(l_3)]+f[tuple(f_3)])*1/280
                )/dx
            )
            
        else:
            raise NotImplementedError("Can only compute 3, 5, 7 or 9 point stencils.")
        
    return gradient

def lapl_discrete(field,dxdydz=(1,1,1)):
    """
    Applies the Laplace operator (\\nabla\\cdot\\nabla) to a field with a standard stencil.
        This function assumes periodic boundary conditions and even spacing
        between points.
    """

    # If dxdydz is not provided, set to tuple of ones.
    try:
        assert len(dxdydz)==field.ndim
    except:
        dxdydz = (1,)*field.ndim

    f = field.copy()

    # Assume periodic boundary conditions by pasting
    # the edge faces on opposite sides of the volume.
    for axis, dimension in enumerate(f.shape):
        
        f_slice, l_slice = (
            [slice(0,d) for d in f.shape],
            [slice(0,d) for d in f.shape],
        )
        
        f_slice[axis] = slice(0,1)
        l_slice[axis] = slice(-1,dimension)
        
        f = np.concatenate((f[tuple(l_slice)], f, f[tuple(f_slice)]), axis=axis)

    # Calculate Laplacian using a standard stencil.
    lapl = 0.
    for axis, (dx, dimension) in enumerate(zip(dxdydz,f.shape)):
        f_slices, l_slices, m_slices = (
            [slice(1,-1) for d in f.shape],
            [slice(1,-1) for d in f.shape],
            [slice(1,-1) for d in f.shape],)
        
        f_slices[axis] = slice(0,-2)
        l_slices[axis] = slice(2,dimension)
        
        lapl += (f[tuple(f_slices)] + f[tuple(l_slices)] - 2*f[tuple(m_slices)]) / dx**2
        
    return lapl

def gradient_FT(field, dxdydz=None):
    """
    Calculates the gradient of a field by Fourier transforming and inverse Fourier transformaing.
    """

    # If dxdydz is not provided, set to tuple of ones.
    try:
        assert len(dxdydz)==field.ndim
    except:
        dxdydz = (1,)*field.ndim
    
    f = field.copy()
    
    kxkykz = 2*np.pi*np.array([np.fft.fftfreq(shape, dx) for (shape,dx) in zip(f.shape,dxdydz)])
    KxKyKz = np.meshgrid(*kxkykz,indexing='ij')
    
    gradient = [ifftn(fftn(f)*Kx*1j).real for Kx in KxKyKz]
    
    return gradient

def invlapl_FT(field, dxdydz=None):
    """
    Given a field f and the equation
            \\nabla^2 V = f
    where \\nabla^2 = \\nabla\\cdot\\nabla is the Laplace operator solve for the potential V.
    """

    # If dxdydz is not provided, set to tuple of ones.
    try:
        assert len(dxdydz)==field.ndim
    except:
        dxdydz = (1,)*field.ndim
    
    # Get frequency space values.
    kxkykz = 2*np.pi*np.array([np.fft.fftfreq(dimension, dx) for (dimension,dx) in zip(field.shape, dxdydz)])
    kxkykz[:,0] = kxkykz[:,1]/1000 # Deal with division by 0
    KxKyKz = np.meshgrid(*kxkykz, indexing='ij')
    K_squared = np.sum([Kx**2 for Kx in KxKyKz], axis=0)
    
    # Calculate Laplacian.
    field_k = fftn(field)
    lapl = -ifftn(field_k/K_squared).real
    
    return lapl

def curl_discrete(vx, vy, vz=None, dxdydz=None):
    """
    Calculates the curl of a vector field. Defined for 2 and 3 dimensions.

    Parameters
    ----------
        vx, vy, vz : np.ndarray
            Field values of a vector field. Must give at least vx and vy.
        dxdydz : iterable
            Step size between different values of field, one iterable for each dimension.

    Returns
    -------
        curl of a vector field.
    """

    # If dxdydz is not provided, set to tuple of ones.
    try:
        assert len(dxdydz)==vx.ndim
    except:
        dxdydz = (1,)*vx.ndim

    # calculate necessary gradients
    dvxdy, dvxdz = gradient_discrete(vx,dxdydz)[1:]
    dvydx, dvydz = gradient_discrete(vy,dxdydz)[::2]
    try:
        dvzdx, dvzdy = gradient_discrete(vz,dxdydz)[:2]
        curlx = dvzdy-dvydz
        curly = dvxdz-dvzdx
        curlz = dvydx-dvxdy
        return curlx, curly, curlz
    except:
        curlz = dvydx-dvxdy
        return curlz