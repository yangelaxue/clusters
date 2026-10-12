"""
Script to animate all raw and derived fields from a PLUTO simulation.

Inputs are:
    Path
    snapshots (optional)
        List of snapshots which to create an animation from.
        Must be 'all' or list of integers.
    slc (optional)
        Which slice to animate.
    los (optional)
        Along which line of sight should a slice be taken.
    c (optional)
        Only plot every c-th data point in the grid.
    L_0 (optional)
        Length unit along the x and y axes.
    t_0 (optional)
        Time unit to label.
TODO Need to plot varticity and curlB
"""

def plot(val,varName,snapshot,vmin=None,vmax=None,norm=None,cbar_label=None):

    # Define x and y values.
    z = get_redshift(snapshot)
    L = IA2.dL * IA2.dim * CGS.pc / (1+z)
    if 'pc' in L_0:
        L /= CGS.pc
    if L_0.startswith('k'):
        L /= 1e3
    elif L_0.startswith('M'):
        L /= 1e6
    else:
        raise NotImplementedError('Does not account for prefixes other than k and M.')
    x, y = (np.linspace(0,L,IA2.dim//c),) * 2

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
    if slc==-1:
        if snapshot not in centres.keys():
            ax.set_title("not through cluster centre",c='r')
    im = ax.pcolormesh(x,y,val,vmin=vmin,vmax=vmax,norm=norm,cmap=cmap)
    cbar = get_cbar(im,fig,ax)
    cbar.ax.set_title(cbar_label,loc='left')

    plt.savefig(os.path.join(savePath,f'{varName}_{snapshot}_slice_los{los}.png'),dpi=DefaultStyle.figkwargs['dpi'])
    plt.close()

def plot_raw():

    for field in IA2.fields:

        vals = ia2.gen_vals(field,snapshots=snapshots,slcs=cents,los=los,c=c)
        
        for s in snapshots:

            # See if the field values loaded.
            try:
                val = next(vals)
            except:
                print(f"... skipping plotting {field} for snapshot {s} ...")
                continue

            print(f"Plotting {field} for snapshot {s}.")
            
            z = get_redshift(s)

            norm = None
            if 'Density' in field:
                val /= CGS.mp
                val *= (1+z)**3
                cbar_label = r'$m_{\rm P}/{\rm cm}^3$'
                norm = colors.LogNorm()
            if field=='Temperature':
                cbar_label = 'K'
            if 'velocity' in field:
                val /= 1e5
                cbar_label=r'km/s'
            if field in ('Bx','By','Bz'):
                val /= 1e-6
                val *= (1+z)**2
                cbar_label = r'$\mu{\rm B}$'

            vmin, vmax = None, None
            if field in ('x-velocity','y-velocity','z-velocity','Bx','By','Bz'):
                vmax = np.max(np.abs(val))
                vmin = -vmax

            plot(val,field,s,vmin,vmax,norm,cbar_label)

def plot_derived():

    # plot v^2 field
    vxs = ia2.gen_vals('x-velocity',snapshots=snapshots,slcs=cents,los=los,c=c)
    vys = ia2.gen_vals('y-velocity',snapshots=snapshots,slcs=cents,los=los,c=c)
    vzs = ia2.gen_vals('z-velocity',snapshots=snapshots,slcs=cents,los=los,c=c)
    for s, in snapshots:
        try:
            vx, vy, vz = next(vxs), next(vys), next(vzs)
        except:
            print(f"... skipping plotting v2 for snapshot {s} ...")
        print(f"Plotting v2 for snapshot {s}.")
        v2 = vx**2 + vy**2 + vz**2
        plot(v2,'v2',s,0,None,None,'$v^2$ cgs units.')        

    # plot ek field
    rhos = ia2.gen_vals('Density',snapshots=snapshots,slcs=cents,los=los,c=c)
    vxs = ia2.gen_vals('x-velocity',snapshots=snapshots,slcs=cents,los=los,c=c)
    vys = ia2.gen_vals('y-velocity',snapshots=snapshots,slcs=cents,los=los,c=c)
    vzs = ia2.gen_vals('z-velocity',snapshots=snapshots,slcs=cents,los=los,c=c)
    for s in snapshots:
        try:
            rho, vx, vy, vz = next(rhos), next(vxs), next(vys), next(vzs)
        except:
            print(f'... skipping plotting v2 for snapshot {s} ...')
        print(f"Plotting ek for snapshot {s}.")
        ek = 0.5 * rho * (vx**2 + vy**2 + vz**2)
        plot(ek,'ek',s,0,None,None,'$v^2$ cgs units.')

    # Plot B2 field
    Bxs = ia2.gen_vals('Bx',snapshots=snapshots,slcs=cents,los=los,c=c)
    Bys = ia2.gen_vals('By',snapshots=snapshots,slcs=cents,los=los,c=c)
    Bzs = ia2.gen_vals('Bz',snapshots=snapshots,slcs=cents,los=los,c=c)
    for s in snapshots:
        try:
            Bx, By, Bz = next(Bxs), next(Bys), next(Bzs)
        except:
            print(f"... skipping plotting B2 for snapshot {s} ...")
        print(f"Plotting B2 for snapshot {s}.")
        B2 = Bx**2 + By**2 + Bz**2
        plot(B2,'B2',s,0,None,None,'$v^2$ cgs units.')

def save_centres(centres):

    snaps, cents = [], []
    for s, centre in centres.items():
        snaps.append([s])
        cents.append(centre)

    arr  = np.concatenate([snaps,cents],axis=1)
    np.savetxt(centrefName,arr,header='snapshot c centre')

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
    parser.add_argument("-slc", help="Which slice to plot.", type=int,required=False)
    parser.add_argument("-los", help="Which line-of-sight to plot a slice from.",required=False)
    parser.add_argument("-c", help="How small a subset to plot.", type=int,required=False)
    parser.add_argument("-L_0", help="Langth unit to plot.", type=int,required=False)
    # parser.add_argument("-t_0", help="Time unit to plot.", type=int,required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    if not args.snapshots:
        snapshots = list(IA2.snapshots)
    else:
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]    
    slc = args.slc if args.slc else -1
    los = args.los if args.los else '0'
    if los in {'0','x','X'}:
        los = 0
    elif los in {'1','y','Y'}:
        los = 1
    elif los in {'2','z','Z'}:
        los = 2
    else:
        raise NotImplementedError("Specified los must be in 0, x, X, 1, y, Y, 2, z, Z")
    c = args.c if args.c else 1
    L_0 = args.L_0 if args.L_0 else 'kpc'

    if 'E3A' in os.getcwd().split('/')[-1] or 'E3A' in Path.split('/')[-1]:
        IA2.dim = 1024

    ia2 = IA2Data(Path)

    savePath = os.path.join(Path,'visualise')
    if not os.path.exists(savePath):
        os.mkdir(savePath)

    # Calculate or load centres
    cents = []
    if slc==-1:
        centres = {}
        centrefName = os.path.join(Path,"centre.txt")
        # Load saved centres
        if os.path.exists(centrefName):
            f = np.loadtxt(centrefName,skiprows=1)
            for s in snapshots:
                for line in f:
                    if s in line:
                        centres.update({s : line[1:]})
                        continue
        # Calculate and save centres
        for s in snapshots:
            if s not in centres:
                try:
                    centre = ia2.get_centre(snapshot=s,c=1)
                    centres.update({s : centre//c})
                    save_centres(centres)
                except:
                    pass
        # Fill in remainder with grid centres
        for s in snapshots:
            if s in centres.keys():
                cents.append(int(centres[s][los]))
            else:
                cents.append(IA2.dim//c//2)
    else:
        for s in snapshots:
            cents.append(slc)

    plot_raw()
    plot_derived()
