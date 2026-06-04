"""
Functions used to help visualise MHD field data.

Classes
-------
    DefaultStyle
        Class that stores my default figure styles parameters.
        If no kwargs are specified, use those defined here.

Functions
---------
    set_style
        Set the figure style to my own.
    get_cbar
        Create and return a colorbar sized to fit an axis.
    get_slice
        When given a 3D volume, return just a slice which to plot onto a meshplot.
    animate_line
        Given a time-ordered list of line data points, plot and save an animated line plot.
    animate_meshgrid
        Given a time-ordered list of mesh points, plot and save an animated mesh plot.

TODO: consider var colors.

Author: Angela
Date: March 2026
"""

import numpy as np
import matplotlib as mpl
from matplotlib import rc, animation
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import cmasher as cmr

# Style class

class DefaultStyle:

    colors = ('b','orange','red','green','orchid','gold')
    markers = ('x','^','D','*','v','+')

    figkwargs = {
        'figsize' : (5,4),
        'dpi' : 300,
        'tight_layout' : True,
    }
    cbarkwargs = {
        'cbar_size' : .2,
        'cbar_pad' : .2,
        'cbar_loc' : 'right',
        'cbar_label' : None,
    }
    animationkwargs = {
        'interval' : 300,
    }

    def get_varkwargs(varName):

        kwargs = {
            'title' : varName,
            'cbar_cmap' : 'viridis',
            'cbar_label' : varName,
            'units' : 'arb. units'
        }

        if varName in {'rho','Density'}:
            kwargs.update({
                'title' : 'Density',
                'cbar_cmap' : 'plasma',
                'cbar_label' : r'$\rho$',
                'units' : r'$m_{\rm p}/{\rm cm}^3$'
            })
        elif varName in {'Dark_Matter_Density',}:
            kwargs.update({
                'title' : 'Dark Matter Density',
                # 'cbar_cmap' : 'plasma', #TODO
                'cbar_label' : r'$\rho$ (DM)',
                'units' : r'$m_{\rm p}/{\rm cm}^3$'
            })
        elif varName in {'prs',}:
            kwargs.update({
                'title' : 'Pressure',
                'cbar_cmap' : 'Spectral_r',
                'cbar_label' : r'$P$',
                # 'units' : r'code units',
            })
        elif varName in ('tem','Temperature'):
            kwargs.update({
                'title' : 'Temperature',
                # 'cbar_cmap' : 'plasma', #TODO
                'cbar_label' : r'$T$',
                'units' : r'K',
            })
        elif varName in {'v',}:
            kwargs.update({
                'title' : 'Velocity',
                'cbar_cmap' : 'Blues',
                'cbar_label' : r'$v$',
                'units' : r'km/s',
            })
        elif varName in {'v2',}:
            kwargs.update({
                'title' : '(twice) Specific Kinetic Energy',
                'cbar_cmap' : 'Blues',
                'cbar_label' : r'$v^2$',
                'units' : r'km$^2$/s$^2$',
            })
        elif varName in {'vor',}:
            kwargs.update({
                'title' : 'Vorticity',
                'cbar_cmap' : 'Blues',
                'cbar_label' : r'$|\nabla\times v|$',
                'units' : '\n' + r'${\rm km}/{\rm s}/{\rm kpc}$', #TODO
            })
        elif 'vx' in varName or 'velocity' in varName:
            kwargs.update({
                'cbar_cmap' : 'bwr',
                'units' : r'km/s',
            })
            if varName in ('vx1','x-velocity'):
                kwargs.update({
                    'cbar_label' : r'$v_x$',
                    'title' : r'$x$-Velocity',
                })
            elif varName in ('vx2','y-velocity'):
                kwargs.update({
                    'cbar_label' : r'$v_y$',
                    'title' : r'$y$-Velocity',
                })
            else:
                kwargs.update({
                    'cbar_label' : r'$v_z$',
                    'title' : r'$z$-Velocity',
                })
        elif varName in {'xraySB','xray'}:
            kwargs.update({
                'title' : 'X-ray Surface Brightness',
                'cbar_cmap' : 'magma',
                'cbar_label' : r'${\rm SB}_{\rm X}$',
                # 'units' : 'arb. units',
            })
        elif varName in {'SZSB','SZ'}:
            kwargs.update({
                'title' : 'SZ Surface Brightness',
                'cbar_cmap' : 'pink',
                'cbar_label' : r'${\rm SB}_{\rm SZ}$',
                # 'units' : 'arb. units',
            })
        elif varName in {'B2',}:
            kwargs.update({
                'title' : '(twice) Magnetic Energy',
                'cbar_cmap' : 'bone', 
                'cbar_label' : r'$B^2$',
                # 'units' : 'arb. units',
            })
        elif varName in {'curlB',}:
            kwargs.update({
                'title' : '$|$Curl of Magnetic Field$|$',
                'cbar_cmap' : 'bone',
                'cbar_label' : r'$|\nabla\times B|$',
                # 'units' : 'arb. units',
            })
        elif varName.startswith('B'):
            kwargs.update({
                'cbar_cmap' : cmr.wildfire,
                # 'units' : 'arb. units',
            })
            if varName in {'Bx1','Bx'}:
                kwargs.update({
                    'title' : r'$x$-Magnetic Field',
                    'cbar_label' : r'$B_x$',
                })
            elif varName in {'Bx2','By'}:
                kwargs.update({
                    'title' : r'$y$-Magnetic Field',
                    'cbar_label' : r'$B_y$',
                })
            elif varName in {'Bx3','Bz'}:
                kwargs.update({
                    'title' : r'$z$-Magnetic Field',
                    'cbar_label' : r'$B_z$',
                })
        else:
            # TODO replace this hardcode with the name of the class.
            print(f'{"DefaultStyle"} catalogue does not include {varName}. Returning to default.')

        return kwargs

# Generally useful functions

def set_style():

    fontsize = 10

    if os.environ['USER']=='yange':
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': fontsize})
        rc('text', usetex=True)
    else:
        rc('font', **{'size': fontsize})
    rc('axes', **{'titlesize': fontsize})

    plt.rcParams['axes.axisbelow'] = True
    mpl.rcParams['lines.linewidth'] = 1

    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
        color=DefaultStyle.colors,
        # marker=DefaultStyle.markers,
        )

def get_cbar(im,fig,ax,**_kwargs):
    """
    Return a cbar of the same size as the plot.
    """

    kwargs = DefaultStyle.cbarkwargs
    kwargs.update(_kwargs)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes(kwargs['cbar_loc'], size=kwargs['cbar_size'], pad=kwargs['cbar_pad'])
    cbar = fig.colorbar(im,cax=cax)
    return cbar

def get_slice(val,slc=0,los='x',c=1):
    """
    From a 3D array, return a 2D array with a specified line-of-sight and slice.
    Can be mesh grids.
    
    Parameters
    ----------
        val : np.ndarray
            Array with 2 or 3 dimensions which to return a slice of.
        slc : int
            Which slice to return of a 3D array.
        los : int, str
            Axis which slice should be taken along.
        c : int
            Return subset of data at regular grid intervals of c.
    """

    if val.ndim==2:
        return val[::c,::c]
    if val.ndim==3:
        if los in [0,'x','X']:
            ret = val[slc,:,:]
        elif los in [1,'y','Y']:
            ret = val[:,slc,:]
        elif los in [2,'z','Z']:
            ret = val[:,:,slc]
        else:
            raise ValueError("Ling-of-sight must be x, y, z, X, Y, Z, 0, 1 or 2.")
        ret = ret[::c,::c]
    else: NotImplementedError("Dimension of val must be 2 or 3.")
    return ret
    
#%% Animations

def animate_line(y_vals,x_vals,fig,ax,labels=[],savePath='lines.mp4',**_kwargs):

    kwargs = DefaultStyle.animationkwargs | DefaultStyle.figkwargs
    kwargs.update(_kwargs)

    # x-points may be static.
    if np.ndim(x_vals)==1:
        x_vals = [x_vals,]*len(y_vals)
    assert len(y_vals)==len(x_vals)
    # Define color
    try:
        color = _kwargs['color']
    except:
        color = None
    # Define marker
    try:
        marker = _kwargs['marker']
    except:
        marker = None

    line, = ax.plot(x_vals[0],y_vals[0],color=color,marker=marker)
    if labels:
        label = ax.annotate(labels[0],(20,20),xycoords='axes pixels',ha='left',va='bottom',color='k',
                            bbox={'facecolor':'white','alpha':0.6,'boxstyle':'square','pad':0.2,'lw':0}
                            )

    y_vals.append(y_vals[-1])
    x_vals.append(x_vals[-1])
    def animate(frame):
        line.set_data(x_vals[frame], y_vals[frame])
        if labels:
            label.set_text(labels[frame])
        return line
    
    anim = animation.FuncAnimation(fig,animate,interval=kwargs['interval'],frames=len(y_vals)-1)
    anim.save(savePath,dpi=kwargs['dpi'])

def animate_meshgrid(vals,fig,ax,varName='',labels=[],savePath='mesh.mp4', **_kwargs):
    """
    Animates a time-ordered list of field value slices.

    Parameters
    ----------
        vals : list of np.ndarray
            Time-ordered list of 2D field values of equal shape.
        varName : str, optional
            Name of variable to plot to retrieve default kwargs.
        labels : list of str, optional
            Labels corresponding to length of vals.
        savePath : str, optional
            Where to save the animation.
        **_kwargs, optional
            x, y : np.ndarray
                Arrays (1D) or meshgrids (2D).
            fig, ax : fig, ax
                Figure and axis which to animate the mesh onto.
    """

    kwargs = DefaultStyle.cbarkwargs | DefaultStyle.animationkwargs | DefaultStyle.figkwargs | DefaultStyle.get_varkwargs(varName)
    kwargs.update(_kwargs)

    # If x and y are not provided... Does not accept X and Y
    try:
        x, y = _kwargs['x'], _kwargs['y']
    except:
        x, y = np.arange(vals[0].shape[0]), np.arange(vals[0].shape[1])
    # Define vmin and vmax
    try:
        vmin, vmax = _kwargs['vmin'], _kwargs['vmax']
        vmin = _kwargs['vmin'] if _kwargs['vmin'] else np.min(vals)
        vmax = _kwargs['vmax'] if _kwargs['vmax'] else np.max(vals)
    except:
        vmin, vmax = np.min(vals), np.max(vals)
        
    im = ax.pcolormesh(x,y,vals[0],vmin=vmin,vmax=vmax,cmap=kwargs['cbar_cmap'])
    if labels:
        label = ax.annotate(labels[0],(20,20),xycoords='axes pixels',ha='left',va='bottom',color='k',
                            bbox={'facecolor':'white','alpha':0.6,'boxstyle':'square','pad':0.2,'lw':0}
                            )
    cbar = get_cbar(im,fig,ax)
    if kwargs['cbar_label']:
        cbar.ax.set_title(kwargs['cbar_label']+' '+kwargs['units'],loc='left')

    vals.append(vals[-1])
    def animate(frame):
        im.set_array(vals[frame].flatten())
        if labels:
            label.set_text(labels[frame])

    anim = animation.FuncAnimation(fig,animate,interval=kwargs['interval'], frames=len(vals)-1)
    anim.save(savePath,dpi=kwargs['dpi'])

    #%% NOT IMPLEMENTED

# def animate_vfield(vx1_list, vx2_list, intv:int=1, **_kwargs):
#     """
#     Animates a 2D field.

#     Parameters
#     ----------
#         vx1_list, vx2_list : list
#             Time ordered list of 2D arrays of the x and y components of a field respectively.
#         kwargs
#             Keyword arguments to specify the appearence of the animation.
#             title : str - title shown on the top of the animation.
#             X, Y : np.ndarray - X and Y coordinates of the same shape of the variables.
#             intv : int - interval between adjacent field lines to be sketched.
#             figsize : tuple - size in inches of the animation.
#             save_dir : str - where and name of the output movie.
#     """
    
#     pp = PlutoPlots
#     pp.kwargs.update(_kwargs)
    
#     if not pp.kwargs['XY']:
#         pp.get_XY(vx1_list[0].shape)
#     X, Y = pp.kwargs['XY']
    
#     fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)
    
#     ax.set_aspect('equal')
#     ax.set_title(f"{pp.kwargs['title']}")
#     ax.set_xlabel(f"{pp.kwargs['xlabel']}")
#     ax.set_ylabel(f"{pp.kwargs['ylabel']}")
    
#     ims = []
#     for i in range(len(vx1_list)):
#         im = ax.quiver(X[::intv,::intv], Y[::intv,::intv], vx1_list[i][::intv,::intv] ,vx2_list[i][::intv,::intv], animated=True)
#         ims.append([im])
    
#     ani = animation.ArtistAnimation(fig, ims, interval=pp.kwargs['interval'])
#     ani.save(pp.kwargs['save_dir'],dpi=150)

# def animate_vstream(vx1_list, vx2_list, density, **_kwargs):

#     # TODO: this doesn't work.

#     """
#     Animates a 2D field.

#     Parameters
#     ----------
#         vx1_list, vx2_list : list
#             Time ordered list of 2D arrays of the x and y components of a field respectively.
#         kwargs
#             Keyword arguments to specify the appearence of the animation.
#             title : str - title shown on the top of the animation.
#             X, Y : np.ndarray - X and Y coordinates of the same shape of the variables.
#             intv : int - interval between adjacent field lines to be sketched.
#             figsize : tuple - size in inches of the animation.
#             save_dir : str - where and name of the output movie.
#     """

#     pp = PlutoPlots
#     pp.kwargs.update(_kwargs)
    
#     if not pp.kwargs['XY']:
#         pp.get_XY(vx1_list[0].shape)
#     X, Y = pp.kwargs['XY']
    
#     fig, ax = plt.subplots(1,1,figsize=pp.kwargs['figsize'],tight_layout=True)
    
#     ax.set_aspect('equal')
#     ax.set_title(f"{pp.kwargs['title']}")
#     ax.set_xlabel(f"{pp.kwargs['xlabel']}")
#     ax.set_ylabel(f"{pp.kwargs['ylabel']}")

#     stream = ax.streamplot(X,Y,vx1_list[0], vx2_list[0], density=density)

#     def animate(i):
#         ax.collections = [] # clear lines streamplot
#         ax.patches = [] # clear arrowheads streamplot
#         # dy = -1 + iter * 0.01 + Y**2
#         # dx = np.ones(dy.shape)
#         # dyu = dy / np.sqrt(dy**2 + dx**2)
#         # dxu = dx / np.sqrt(dy**2 + dx**2)
#         stream = ax.streamplot(X,Y,vx1_list[i], vx2_list[i], density=2,arrowsize=1)
#         print(i)
#         return stream
    
#     # ims = []
#     # for i in range(len(vx1_list)):
#     #     im = ax.streamplot(X, Y, vx1_list[i] ,vx2_list[i])
#     #     ims.append([im])
    
#     ani = animation.FuncAnimation(fig, animate, interval=pp.kwargs['interval'])
#     ani.save(pp.kwargs['save_dir'],dpi=150)
