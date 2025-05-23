"""
Functions used to help visualise data.
In theory, this should be designed to work with any type of data.

Author: Angela Xue
"""

#%% Imports

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

#################
#%% Default Class
#################

class PlutoPlots:
    """
    Define default matplitlib settings for PLUTO plots.
    """

    kwargs = {
        # Axes
        'XY' : None,
        # Labels
        'title' : None,
        'xlabel' : None,
        'ylabel' : None,
        # Cbar
        'vmin' : None,
        'vmax' : None,
        'cbar:cmap' : 'plasma',
        'cbar:label' : None,
        'cbar:size' : .2,
        'cbar:pad' : .2,
        'cbar:loc' : 'right',
        # Plots
        'figsize' : (6,6),
        'save_dir' : './movie.mp4',
        'interval' : 300,
    }

    def set_style(self,):

        from matplotlib import rc
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
        rc('text', usetex=True)
        rc('axes', **{'titlesize': 10})

        plt.rcParams['axes.axisbelow'] = True

    def get_XY(self,shape):
        """
        Generate a default meshgrid.
        """
        X, Y = np.meshgrid(np.linspace(0,1.,shape[0]), np.linspace(0,1.,shape[1]), indexing='xy')
        self.kwargs['XY'] = X, Y


def get_cbar(im,fig,ax,**_kwargs):
    """
    Return a cbar of the same size as the plot.
    """

    pp = PlutoPlots()
    pp.kwargs.update(_kwargs)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes(pp.kwargs['cbar:loc'], size=pp.kwargs['cbar:size'], pad=pp.kwargs['cbar:pad'])
    cbar = fig.colorbar(im,cax=cax)
    return cbar
    
def get_fig_ax(**_kwargs):

    pp = PlutoPlots()
    pp.kwargs.update(_kwargs)

    fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)
    
    ax.set_aspect('equal')
    ax.set_xlabel(pp.kwargs['x_label'])
    ax.set_ylabel(pp.kwargs['y_label'])


    return fig, ax

####################
#%% Helper Functions
####################

def get_list2D(var_list,s=0):
    """
    From list of 3D arrays, return a list of 2D arrays from a specified slice s.
    """
    return [_var[s] for _var in var_list] if var_list[0].ndim==3 else var_list

# def set_plot(**_kwargs):

#     pp = PlutoPlots()
#     kwargs = pp.kwargs
#     kwargs.update(_kwargs)

#     fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)


############
#%% Animate!
############

def animate_vars(var_list, **_kwargs):
    """
    Animates a 2D grid of values.

    Parameters
    ----------
        var_list : list
            Time ordered list of 2D arrays which are to be plotted.
        _kwargs
            Keyword arguments to specify the appearence of the animation.
            See PlutoPlots.
        """
    
    pp = PlutoPlots()
    pp.kwargs.update(_kwargs)
    
    if not pp.kwargs['XY']:
        pp.get_XY(var_list[0].shape)
    X, Y = pp.kwargs['XY']
    
    fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)

    ax.set_aspect('equal')
    
    ims = []
    for i in range(len(var_list)):
        im = ax.pcolormesh(X, Y, var_list[i], animated=True, cmap=pp.kwargs['cmap'], vmin=pp.kwargs['vmin'], vmax=pp.kwargs['vmax'])
        ims.append([im])
    
    if pp.kwargs['vmin']!=None and pp.kwargs['vmax']!=None:
        cbar = get_cbar(im,fig,ax)
        if pp.kwargs['cbar_label']:
            cbar.set_label(pp.kwargs['cbar_label'])
    
    ani = animation.ArtistAnimation(fig, ims, interval=pp.kwargs['interval'])
    ani.save(pp.kwargs['save_dir'],dpi=150)

def animate_vfield(vx1_list, vx2_list, intv:int=1, **_kwargs):
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
    
    pp = PlutoPlots
    pp.kwargs.update(_kwargs)
    
    if not pp.kwargs['XY']:
        pp.get_XY(vx1_list[0].shape)
    X, Y = pp.kwargs['XY']
    
    fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)
    
    ax.set_aspect('equal')
    ax.set_title(f"{pp.kwargs['title']}")
    ax.set_xlabel(f"{pp.kwargs['xlabel']}")
    ax.set_ylabel(f"{pp.kwargs['ylabel']}")
    
    ims = []
    for i in range(len(vx1_list)):
        im = ax.quiver(X[::intv,::intv], Y[::intv,::intv], vx1_list[i][::intv,::intv] ,vx2_list[i][::intv,::intv], animated=True)
        ims.append([im])
    
    ani = animation.ArtistAnimation(fig, ims, interval=pp.kwargs['interval'])
    ani.save(pp.kwargs['save_dir'],dpi=150)

def animate_vstream(vx1_list, vx2_list, density, **_kwargs):

    # TODO: this doesn't work.

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

    pp = PlutoPlots
    pp.kwargs.update(_kwargs)
    
    if not pp.kwargs['XY']:
        pp.get_XY(vx1_list[0].shape)
    X, Y = pp.kwargs['XY']
    
    fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)
    
    ax.set_aspect('equal')
    ax.set_title(f"{pp.kwargs['title']}")
    ax.set_xlabel(f"{pp.kwargs['xlabel']}")
    ax.set_ylabel(f"{pp.kwargs['ylabel']}")

    stream = ax.streamplot(X,Y,vx1_list[0], vx2_list[0], density=density)

    def animate(i):
        ax.collections = [] # clear lines streamplot
        ax.patches = [] # clear arrowheads streamplot
        # dy = -1 + iter * 0.01 + Y**2
        # dx = np.ones(dy.shape)
        # dyu = dy / np.sqrt(dy**2 + dx**2)
        # dxu = dx / np.sqrt(dy**2 + dx**2)
        stream = ax.streamplot(X,Y,vx1_list[i], vx2_list[i], density=2,arrowsize=1)
        print(i)
        return stream
    
    # ims = []
    # for i in range(len(vx1_list)):
    #     im = ax.streamplot(X, Y, vx1_list[i] ,vx2_list[i])
    #     ims.append([im])
    
    ani = animation.FuncAnimation(fig, animate, interval=pp.kwargs['interval'])
    ani.save(pp.kwargs['save_dir'],dpi=150)