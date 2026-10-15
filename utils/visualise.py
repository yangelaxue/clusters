"""
Functions used to help visualise MHD field data.

Classes
-------
    Style
        Class that stores my default figure styles parameters.
        Methods
        -------
            set
                sets style of plots.
            get_cbar
                Returns a colorbar for a plot.
            get_varkwargs
                Returns variable specific kwargs.

Functions
---------
    get_slice
        When given a 3D volume, return just a slice which to plot onto a meshplot.

Animations
----------
    # animate_line #TODO
    #     Given a time-ordered list of line data points, plot and save an animated line plot.
    # animate_meshgrid #TODO
    #     Given a time-ordered list of mesh points, plot and save an animated mesh plot.

Author: Angela
Date: March 2026
Revised: October 2026
"""

class Style:

    figkwargs = {
        "dpi" : 300,
        "tight_layout" : True,
    }

    def set():
        import socket
        import matplotlib as mpl

        fontsize = 10
        if socket.gethostname()=='yange':
            mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
            mpl.rc('text', usetex=True)
        mpl.rc('axes', **{'titlesize': fontsize})
        mpl.rc('font', **{'size': fontsize})

        mpl.pyplot.rcParams['axes.axisbelow'] = True
        mpl.rcParams['lines.linewidth'] = 1

        mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
            color=('b','orange','red','green','orchid','gold'),
            )
        
    def get_cbar(im,fig,ax,**_kwargs):
        """ Return a cbar of the same size as the plot. """
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        kwargs = {
            "size" : .2,
            "pad" : 0.,
            "loc" : "right",
            }
        kwargs.update(_kwargs)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes(kwargs["loc"], size=kwargs["size"], pad=kwargs["pad"])
        cbar = fig.colorbar(im,cax=cax)
        return cbar

    def get_varkwargs(varName):

        # Default kwargs for unidentified varNames.
        kwargs = {
            "title" : varName,
            "cmap" : "viridis",
            "symbol" : varName
        }

        if varName in {"rho","Density"}:
            kwargs.update({
                "title" : "Density",
                "cmap" : "plasma",
                "symbol" : r"$\rho$",
            })
        elif varName in {"dm","Dark_Matter_Density"}:
            kwargs.update({
                "title" : "Dark Matter Density",
                # "cmap" : "plasma", #TODO
                "symbol" : r"$\rho_{\rm DM}$",
            })
        elif varName in {"prs"}:
            kwargs.update({
                "title" : "Pressure",
                "cmap" : "Spectral_r",
                "symbol" : "$P$",
            })
        elif varName in {"tem","Temperature"}:
            kwargs.update({
                "title" : "Temperature",
                # "cmap" : "plasma", #TODO
                "symbol" : "T$",
            })
        elif varName in {"vx1","vx","x-velocity"}:
            kwargs.update({
                "title" : "$x$-Velocity",
                "cmap" : "bwr",
                "symbol" : "$v_x$",
            })
        elif varName in {"vx2","vy","y-velocity"}:
            kwargs.update({
                "title" : "$y$-Velocity",
                "cmap" : "bwr",
                "symbol" : "$v_y$",
            })
        elif varName in {"vx3","vz","z-velocity"}:
            kwargs.update({
                "title" : "$z$-Velocity",
                "cmap" : "bwr",
                "symbol" : "$v_z$",
            })
        elif varName in {"v2/2",}:
            kwargs.update({
                "title" : "Specific Kinetic Energy",
                # 'cbar_cmap' : 'Blues', #TODO
                "symbol" : "$v^2/2$",
            })
        elif varName.lower() in {"KE","kinetic"}:
            kwargs.update({
                "title" : "Kinetic Energy",
                # 'cbar_cmap' : 'Blues', #TODO
                "symbol" : "$E_K$",
            })
        elif varName in {"vor","vorticity"}:
            kwargs.update({
                "title" : "Vorticity",
                # "cbar_cmap" : "Blues", #TODO
                "symbol" : r"$|\nabla\times v|$",
            })
        elif varName in {"Bx1","Bx"}:
            kwargs.update({
                "title" : "$x$-Magnetic Field",
                "cmap" : "BrBG",
                "symbol" : "$B_x$",
            })
        elif varName in {"Bx2","By"}:
            kwargs.update({
                "title" : "$y$-Magnetic Field",
                "cmap" : "BrBG",
                "symbol" : "$B_y$",
            })
        elif varName in {"Bx3","Bz"}:
            kwargs.update({
                "title" : "$z$-Magnetic Field",
                "cmap" : "BrBG",
                "symbol" : "$B_z$",
            })
        elif varName in {"b2/2",}:
            kwargs.update({
                "title" : "Magnetic Energy",
                # 'cbar_cmap' : 'Blues', #TODO
                "symbol" : "$B^2/2$",
            })
        elif varName in {"curlB"}:
            kwargs.update({
                "title" : "Curl Mangetic Field",
                # "cbar_cmap" : "Blues", #TODO
                "symbol" : r"$|\nabla\times B|$",
            })
        elif varName in {"xraySB","xray"}:
            kwargs.update({
                "title" : "X-ray Surface Brightness",
                "cmap" : "bone", #TODO
                "symbol" : r"${\rm SB}_{\rm SZ}$",
            })
        elif varName in {"SZSB","SZ"}:
            kwargs.update({
                "title" : "SZ Surface Brightness",
                "cmap" : "pink", #TODO
                "symbol" : r"${\rm SB}_{\rm SZ}$",
            })
        return kwargs

def get_slice(val,slc=0,los='x',c=1):
    """
    From a 3D array, return a 2D array with a specified line-of-sight and slice.
    """

    if val.ndim==2:
        ret = val[::c,::c]
    if val.ndim==3:
        if los in [0,'x','X']:
            ret = val[slc,:,:]
        elif los in [1,'y','Y']:
            ret = val[:,slc,:]
        elif los in [2,'z','Z']:
            ret = val[:,:,slc]
        else:
            raise ValueError("Line-of-sight must be x, y, z, X, Y, Z, 0, 1 or 2.")
        ret = ret[::c,::c]
    else: NotImplementedError("Dimension of val must be 2 or 3.")
    return ret

#%% Animations!

# def animate_line(y_vals,x_vals,fig,ax,labels=[],savePath='lines.mp4',**_kwargs):

#     kwargs = DefaultStyle.animationkwargs | DefaultStyle.figkwargs
#     kwargs.update(_kwargs)

#     # x-points may be static.
#     if np.ndim(x_vals)==1:
#         x_vals = [x_vals,]*len(y_vals)
#     assert len(y_vals)==len(x_vals)
#     # Define color
#     try:
#         color = _kwargs['color']
#     except:
#         color = None
#     # Define marker
#     try:
#         marker = _kwargs['marker']
#     except:
#         marker = None

#     line, = ax.plot(x_vals[0],y_vals[0],color=color,marker=marker)
#     if labels:
#         label = ax.annotate(labels[0],(20,20),xycoords='axes pixels',ha='left',va='bottom',color='k',
#                             bbox={'facecolor':'white','alpha':0.6,'boxstyle':'square','pad':0.2,'lw':0}
#                             )

#     y_vals.append(y_vals[-1])
#     x_vals.append(x_vals[-1])
#     def animate(frame):
#         line.set_data(x_vals[frame], y_vals[frame])
#         if labels:
#             label.set_text(labels[frame])
#         return line
    
#     anim = animation.FuncAnimation(fig,animate,interval=kwargs['interval'],frames=len(y_vals)-1)
#     anim.save(savePath,dpi=kwargs['dpi'])

# def animate_meshgrid(vals,fig,ax,varName='',labels=[],savePath='mesh.mp4', **_kwargs):
#     """
#     Animates a time-ordered list of field value slices.

#     Parameters
#     ----------
#         vals : list of np.ndarray
#             Time-ordered list of 2D field values of equal shape.
#         varName : str, optional
#             Name of variable to plot to retrieve default kwargs.
#         labels : list of str, optional
#             Labels corresponding to length of vals.
#         savePath : str, optional
#             Where to save the animation.
#         **_kwargs, optional
#             x, y : np.ndarray
#                 Arrays (1D) or meshgrids (2D).
#             fig, ax : fig, ax
#                 Figure and axis which to animate the mesh onto.
#     """

#     kwargs = DefaultStyle.cbarkwargs | DefaultStyle.animationkwargs | DefaultStyle.figkwargs | DefaultStyle.get_varkwargs(varName)
#     kwargs.update(_kwargs)

#     # If x and y are not provided... Does not accept X and Y
#     try:
#         x, y = _kwargs['x'], _kwargs['y']
#     except:
#         x, y = np.arange(vals[0].shape[0]), np.arange(vals[0].shape[1])
#     # Define vmin and vmax
#     try:
#         vmin, vmax = _kwargs['vmin'], _kwargs['vmax']
#         vmin = _kwargs['vmin'] if _kwargs['vmin'] else np.min(vals)
#         vmax = _kwargs['vmax'] if _kwargs['vmax'] else np.max(vals)
#     except:
#         vmin, vmax = np.min(vals), np.max(vals)
        
#     im = ax.pcolormesh(x,y,vals[0],vmin=vmin,vmax=vmax,cmap=kwargs['cbar_cmap'])
#     if labels:
#         label = ax.annotate(labels[0],(20,20),xycoords='axes pixels',ha='left',va='bottom',color='k',
#                             bbox={'facecolor':'white','alpha':0.6,'boxstyle':'square','pad':0.2,'lw':0}
#                             )
#     cbar = get_cbar(im,fig,ax)
#     if kwargs['cbar_label']:
#         cbar.ax.set_title(kwargs['cbar_label']+' '+kwargs['units'],loc='left')

#     vals.append(vals[-1])
#     def animate(frame):
#         im.set_array(vals[frame].flatten())
#         if labels:
#             label.set_text(labels[frame])

#     anim = animation.FuncAnimation(fig,animate,interval=kwargs['interval'], frames=len(vals)-1)
#     anim.save(savePath,dpi=kwargs['dpi'])

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
