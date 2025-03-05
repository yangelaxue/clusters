"""
Functions used to help visualise data.
In theory, this should be designed to work with any type of data.

Author: Angela Xue
"""

#%% Imports

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

#%% Begin.

def animate_vars(var_list, **_kwargs):
    """
    Animates a 2D grid of values.

    Parameters
    ----------
        var_list : list
            Time ordered list of 2D arrays which are to be plotted.
        kwargs
            Keyword arguments to specify the appearence of the animation.
            title : str - title shown on the top of the animation.
            XY : tuple - X and Y coordinates of the same shape of the variables.
            vmin, vmax : float - respectively, the minimum and maximum numbers to be shown throughout the animation.
            figsize : tuple - size in inches of the animation.
            cmap : str - colors used to plot the variables.
            save_dir : str - where and name of the output movie.
    """

    kwargs = {
        'title' : None,
        'XY' : None,
        'vmin' : None,
        'vmax' : None,
        'figsize' : (6,6),
        'cmap' : 'plasma',
        'save_dir' : "./movie.mp4",
    }
    kwargs.update(_kwargs)

    title = kwargs['title']
    XY = kwargs['XY']
    if XY:
        X, Y = XY
    else:
        X, Y = np.meshgrid(np.arange(0,var_list[0].shape[0]), np.arange(0,var_list[0].shape[1]), indexing='xy')
    vmin, vmax = kwargs['vmin'], kwargs['vmax']
    cmap = kwargs['cmap']
    figsize = kwargs['figsize']
    save_dir = kwargs['save_dir']
    
    fig, ax = plt.subplots(1,1,figsize=figsize,tight_layout=True)
    
    ax.set_title(f"{title}")
    ax.set_aspect('equal')
    
    ims = []
    for i in range(len(var_list)):
        mesh = ax.pcolormesh(X, Y, var_list[i], animated=True, cmap=cmap, vmin=vmin, vmax=vmax)
        ims.append([mesh])
    
    if vmin!=None and vmax!=None:
        fig.colorbar(mesh,shrink=.8)
    
    ani = animation.ArtistAnimation(fig, ims, interval=300)
    ani.save(save_dir,dpi=150)

def animate_vfield(vx1_list, vx2_list, **_kwargs):
    """
    Animates a 2D field.

    Parameters
    ----------
        vx1_list, vx2_list : list
            Time ordered list of 2D arrays of the x and y components of a field respectively.
        kwargs
            Keyword arguments to specify the appearence of the animation.
            title : str - title shown on the top of the animation.
            X, Y : np.ndarray - X and Y coordinates of the same shape of the variables.
            intv : int - interval between adjacent field lines to be sketched.
            figsize : tuple - size in inches of the animation.
            save_dir : str - where and name of the output movie.
    """

    kwargs = {
        'title' : None,
        'XY' : None,
        'intv' : 4,
        'figsize' : (6,6),
        'save_dir' : "./movie.mp4",
    }
    kwargs.update(_kwargs)

    title = kwargs['title']
    XY = kwargs['XY']
    if XY:
        X, Y = XY
    else:
        X, Y = np.meshgrid(np.arange(0,vx1_list[0].shape[0]), np.arange(0,vx1_list[0].shape[1]), indexing='xy')
    intv = kwargs['intv']
    figsize = kwargs['figsize']
    save_dir = kwargs['save_dir']
    
    fig, ax = plt.subplots(1,1,figsize=figsize,tight_layout=True)
    
    ax.set_title(f"{title}")
    ax.set_aspect('equal')
    
    ims = []
    for i in range(len(vx1_list)):
        im = ax.quiver(X[::intv,::intv], Y[::intv,::intv], vx1_list[i][::intv,::intv] ,vx2_list[i][::intv,::intv], animated=True)
        ims.append([im])
    
    ani = animation.ArtistAnimation(fig, ims, interval=300)
    ani.save(save_dir,dpi=150)