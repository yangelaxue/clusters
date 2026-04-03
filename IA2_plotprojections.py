"""
Script to plot all raw and derived fields from a PLUTO simulation.

Inputs are:
    Path
    snapshots (optional)
        List of snapshots which to create an animation from.
        Must be 'all' or list of integers.
    los (optional)
        Along which line of sight should a projection be taken.
    L_0 (optional)
        Length unit along the x and y axes.
    t_0 (optional)
        Time unit to label.
TODO Need to plot varticity and curlB
"""

from scipy.integrate import trapezoid

def plot(val,varName,snapshot,vmin=None,vmax=None,norm=None,cbar_label=None):

    # Define x and y values.
    z = get_redshift(snapshot)
    L = IA2.L * CGS.pc / (1+z)
    if 'pc' in L_0:
        L /= CGS.pc
    if L_0.startswith('k'):
        L /= 1e3
    elif L_0.startswith('M'):
        L /= 1e6
    else:
        raise NotImplementedError('Does not account for prefixes other than k and M.')
    x, y = (np.linspace(0,L,IA2.dim),) * 2

    # Start plotting.
    cmap = DefaultStyle.get_varkwargs(varName)['cbar_cmap']

    fig, ax = plt.subplots(**DefaultStyle.figkwargs)
    ax.set_aspect('equal')
    ax.set_xlabel(L_0)
    ax.set_ylabel(L_0)
    ax.annotate(
        f'z={z}',(20,20),xycoords='axes pixels',ha='left',va='bottom',color='k',
        bbox={'facecolor':'white','alpha':0.6,'boxstyle':'square','pad':0.2,'lw':0}
        )
    im = ax.pcolormesh(x,y,val,vmin=vmin,vmax=vmax,norm=norm,cmap=cmap)
    cbar = get_cbar(im,fig,ax)
    cbar.ax.set_title(cbar_label,loc='left')

    plt.savefig(os.path.join(savePath,f'{varName}_{snapshot}_projection_los{los}.png'),dpi=DefaultStyle.figkwargs['dpi'])
    plt.close()

def plot_raw():

    for field in IA2.fields:

        vals = ia2.gen_vals(field,snapshots=snapshots)

        for s, val in zip(snapshots,vals):
            
            # See if the field values loaded.
            try:
                assert type(val)==np.ndarray
            except:
                print(f"... skipping plotting {field} for snapshot {s} ...")
                continue
            
            z = get_redshift(s)

            norm = None
            if 'Density' in field:
                val /= CGS.mp
                val *= (1+z)**3
                cbar_label = r'$m_{\rm P}/{\rm cm}^3\times{\rm cm}$'
                norm = colors.LogNorm()
            elif field=='Temperature':
                cbar_label = r'K$\times{\rm cm}$'
            elif 'velocity' in field:
                val /= 1e5
                cbar_label=r'km/s$\times{\rm cm}$'
            elif field in ('Bx','By','Bz'):
                val /= 1e-6
                val *= (1+z)**2
                cbar_label = r'$\mu{\rm B}\times{\rm cm}$'
            else:
                print(f"... skipping plotting {field} for snapshot {s} ...")
                continue

            print(f"Plotting {field} for snapshot {s}.")

            z = get_redshift(s)
            L = IA2.L * CGS.pc / (1+z)
            if 'pc' in L_0:
                L /= CGS.pc
            if L_0.startswith('k'):
                L /= 1e3
            elif L_0.startswith('M'):
                L /= 1e6
            else:
                raise NotImplementedError('Does not account for prefixes other than k and M.')

            val = trapezoid(val,dx=L,axis=los)
            
            vmin, vmax = None, None
            if field in ('x-velocity','y-velocity','z-velocity','Bx','By','Bz'):
                vmax = np.max(np.abs(val))
                vmin = -vmax

            plot(val,field,s,vmin,vmax,norm,cbar_label)



# def plot_raw(snapshot):
#     """ Animate all the raw fields. """

#     for field in IA2.fields:

#         norm = None

#         try: # Check if we have the file downloaded.
#             val = ia2.get_val(field,snapshot=snapshot,redshift=z)
#         except:
#             print(f'... skipping plotting {field} for snapshot {snapshot} ...')
#             continue

#         print(f"Plotting {field} for snapshot {snapshot}.")
#         # Plot everything in predetermined units.
#         if 'Density' in field:
#             val /= CGS.mp
#             val *= (1+z)**3
#             cbar_label = r'$m_{\rm P}/{\rm cm}^3\times{\rm cm}$'
#             norm = colors.LogNorm()
#         elif field=='Temperature':
#             cbar_label = r'K$\times{\rm cm}$'
#         # elif 'velocity' in field:
#         #     val /= 1e5
#         #     cbar_label=r'km/s$\times{\rm cm}$'
#         # elif field in ('Bx','By','Bz'):
#         #     val /= 1e-6
#         #     val *= (1+z)**2
#         #     cbar_label = r'$\mu{\rm B}\times{\rm cm}$'
#         else:
#             continue

#         vmin, vmax = None, None
#         if field in ('x-velocity','y-velocity','z-velocity','Bx','By','Bz'):
#             vmax = np.max(np.abs(val))
#             vmin = -vmax

#         val = trapezoid(val,dx=L/IA2.dim*CGS.pc,axis=los)

#         plot(val,field,snapshot,vmin,vmax,norm,cbar_label)

def plot_derived(snapshot):

    # Plot v^2 field
    try: # Check if we have the file downloaded.
        vx = ia2.get_val('x-velocity',snapshot=snapshot,redshift=z)
        vy = ia2.get_val('y-velocity',snapshot=snapshot,redshift=z)
        vz = ia2.get_val('z-velocity',snapshot=snapshot,redshift=z)
        print(f"Plotting v2 for snapshot = {snapshot}")
        v2 = vx**2 + vy**2 + vz**2
        val = trapezoid(v2,dx=L/IA2.dim*CGS.pc,axis=los)
        plot(val,'v2',snapshot,0,None,None,'$v^2\times{\rm cm}$ cgs units.')
    except:
        print(f'... skipping plotting v2 for snapshot {snapshot} ...')
    try:
        rho = ia2.get_val('Density',snapshot=snapshot,redshift=z,slc=centre[los],los=los,c=c)
        ek = 0.5 * rho * (vx**2 + vy**2 + vz**2)
        print(f"Plotting ek for snapshot = {snapshot}")
        plot(ek,'ek',snapshot,0,None,None,'KE cgs units.')
    except:
        print(f'... skipping plotting ek for snapshot {snapshot} ...')

    # Plot B^2 field
    try:
        Bx = ia2.get_val('Bx',snapshot=snapshot,redshift=z)
        By = ia2.get_val('By',snapshot=snapshot,redshift=z)
        Bz = ia2.get_val('Bz',snapshot=snapshot,redshift=z)
        print(f"Plotting B2 for snapshot = {snapshot}")
        B2 = Bx**2 + By**2 + Bz**2
        val = trapezoid(B2,dx=L/IA2.dim*CGS.pc,axis=los)
        plot(val,'B2',snapshot,0,None,None,'$B^2$ cgs units.')
    except:
        print(f'... skipping plotting B2 for snapshot {snapshot} ...')

if __name__=="__main__":

    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import colors
    import os, argparse

    from utils.enzo import IA2Data, IA2, get_redshift
    from utils.units import CGS
    from utils.visualize import DefaultStyle, get_cbar, set_style

    set_style()

    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("-snapshots", help="Which snapshots to animate.",required=False)
    parser.add_argument("-los", help="Which line-of-sight to plot a projection from.",required=False)
    parser.add_argument("-L_0", help="Langth unit to plot.", type=int,required=False)
    # parser.add_argument("-t_0", help="Time unit to plot.", type=int,required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    if not args.snapshots:
        snapshots = IA2.snapshots
    else:
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]    
    los = args.los if args.los else '0'
    if los in {'0','x','X'}:
        los = 0
    elif los in {'1','y','Y'}:
        los = 1
    elif los in {'2','z','Z'}:
        los = 2
    else:
        raise NotImplementedError("Specified los must be in 0, x, X, 1, y, Y, 2, z, Z")
    L_0 = args.L_0 if args.L_0 else 'kpc'

    if 'E3A' in os.getcwd().split('/')[-1] or 'E3A' in Path.split('/')[-1]:
        IA2.dim = 1024

    ia2 = IA2Data(Path)

    savePath = os.path.join(Path,'visualise')
    if not os.path.exists(savePath):
        os.mkdir(savePath)

    with open(os.path.join(Path,"centre.txt"),"a") as f:
        f.write("snapshot c centre\n")

    for s in snapshots:

        

        plot_raw(s)
        plot_derived(s)
