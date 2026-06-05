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
"""

def plot_raw():
    """ Animate all the raw fields. """

    for var in pd.varNames:

        vals = pd.get_vals(var,snapshots,slc=slc,los=los,c=c,units=units)

        vmin, vmax = None, None
        if 'vx' in var or 'Bx' in var:
            vmax = np.max(np.abs(vals))
            vmin = -vmax

        fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'])
        ax.set_aspect('equal')
        ax.set_xlabel('$x$' + ' ' + L_0)
        ax.set_ylabel('$y$' + ' ' + L_0)

        animate_meshgrid(vals,fig,ax,var,labels,savePath=os.path.join(savePath,var+'_slice.mp4'),x=x,y=y,vmin=vmin,vmax=vmax,units=units)
    
def plot_derived():
    """
    Animtate derived fields
        v^2, vorticity, B^2, curl(B).
    """
    # Animate derived v2 field
    vx1s = pd.gen_vals('vx1',snapshots,slc=slc,los=los,c=c,units=units)
    vx2s = pd.gen_vals('vx2',snapshots,slc=slc,los=los,c=c,units=units)
    if pd.ndim==2:
        vals = [vx1**2 + vx2**2 for vx1,vx2 in zip(vx1s,vx2s)]
    elif pd.ndim==3:
        vx3s = pd.gen_vals('vx3',snapshots,slc=slc,los=los,c=c,units=units)
        vals = [vx1**2 + vx2**2+vx3**2 for vx1,vx2,vx3 in zip(vx1s,vx2s,vx3s)]
    fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'])
    ax.set_aspect('equal')
    ax.set_xlabel('$x$' + ' ' + L_0)
    ax.set_ylabel('$y$' + ' ' + L_0)
    animate_meshgrid(vals,fig,ax,'v2',labels,savePath=os.path.join(savePath,'v2'+'_slice.mp4'),x=x,y=y,vmin=0,vmax=np.max(vals),units=units)

    # Animate derived vor v field
    dxdydz = pd.get_dxdydz(L_0)
    dxdydz = tuple(dx*c for dx in dxdydz)
    vx1s = pd.gen_vals('vx1',snapshots,c=c,units=units)
    vx2s = pd.gen_vals('vx2',snapshots,c=c,units=units)
    if pd.ndim==2:
        vor = (curl_discrete(vx1,vx2,dxdydz=dxdydz) for vx1,vx2 in zip(vx1s,vx2s))
        vormag = ((vorx**2 + vory**2)**.5 for vorx,vory in vor)
    elif pd.ndim==3:
        vx3s = pd.gen_vals('vx3',snapshots,c=c,units=units)
        vor = [curl_discrete(vx1,vx2,vx3,dxdydz=dxdydz) for vx1,vx2,vx3 in zip(vx1s, vx2s, vx3s)]
        vormag = ((vorx**2 + vory**2 + vorz**2)**.5 for vorx,vory,vorz in vor)
    vals = [get_slice(val,slc,los) for val in vormag]
    fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'])
    ax.set_aspect('equal')
    ax.set_xlabel('$x$' + ' ' + L_0)
    ax.set_ylabel('$y$' + ' ' + L_0)
    animate_meshgrid(vals,fig,ax,'vor',labels,savePath=os.path.join(savePath,'vor'+'_slice.mp4'),x=x,y=y,vmin=0,vmax=np.max(vals),units=units)

    if pd.PHYSICS=='MHD':
        # Animate derived B2 field
        Bx1s = pd.gen_vals('Bx1',snapshots,slc=slc,los=los,c=c,units=units)
        Bx2s = pd.gen_vals('Bx2',snapshots,slc=slc,los=los,c=c,units=units)
        if pd.ndim==2:
            vals = [Bx1**2 + Bx2**2 for Bx1,Bx2 in zip(Bx1s,Bx2s)]
        elif pd.ndim==3:
            Bx3s = pd.gen_vals('Bx3',snapshots,slc=slc,los=los,c=c,units=units)
            vals = [Bx1**2 + Bx2**2 + Bx3**2 for Bx1,Bx2,Bx3 in zip(Bx1s,Bx2s,Bx3s)]
        fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'])
        ax.set_aspect('equal')
        ax.set_xlabel('$x$' + ' ' + L_0)
        ax.set_ylabel('$y$' + ' ' + L_0)
        animate_meshgrid(vals,fig,ax,'B2',labels,savePath=os.path.join(savePath,'B2'+'_slice.mp4'),x=x,y=y,vmin=0,vmax=np.max(vals),units=units)

        # Animate derived curlB field
        dxdydz = pd.get_dxdydz(L_0)
        dxdydz = tuple(dx*c for dx in dxdydz)
        Bx1s = pd.gen_vals('Bx1',snapshots,c=c,units=units)
        Bx2s = pd.gen_vals('Bx2',snapshots,c=c,units=units)
        if pd.ndim==2:
            vor = (curl_discrete(Bx1,Bx2,dxdydz=dxdydz) for Bx1,Bx2 in zip(Bx1s,Bx2s))
            vormag = ((vorx**2 + vory**2)**.5 for vorx,vory in vor)
        elif pd.ndim==3:
            Bx3s = pd.gen_vals('Bx3',snapshots,c=c,units=units)
            vor = [curl_discrete(Bx1,Bx2,Bx3,dxdydz=dxdydz) for Bx1,Bx2,Bx3 in zip(Bx1s, Bx2s, Bx3s)]
            vormag = ((vorx**2 + vory**2 + vorz**2)**.5 for vorx,vory,vorz in vor)
        vals = [get_slice(val,slc,los) for val in vormag]
        fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'])
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ kpc')
        ax.set_ylabel('$y$ kpc')
        animate_meshgrid(vals,fig,ax,'curlB',labels,savePath=os.path.join(savePath,'curlB'+'_slice.mp4'),x=x,y=y,vmin=0,vmax=np.max(vals),units=units)

if __name__=="__main__":

    import numpy as np
    from matplotlib import pyplot as plt
    import argparse, os

    from utils.pluto import PlutoData
    from utils.gradient import curl_discrete
    from utils.visualize import set_style, animate_meshgrid, DefaultStyle, get_slice

    set_style()

    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("-snapshots", help="Which snapshots to animate.",required=False)
    parser.add_argument("-slc", help="Which slice to plot.", type=int,required=False)
    parser.add_argument("-units", help="Which units to plot in.", type=str,required=False)
    parser.add_argument("-los", help="Which line-of-sight to plot a slice from.",required=False)
    parser.add_argument("-c", help="How small a subset to plot.", type=int,required=False)
    parser.add_argument("-L_0", help="Langth unit to plot.", type=str,required=False)
    parser.add_argument("-t_0", help="Time unit to plot.", type=str,required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    if not args.snapshots:
        snapshots = 'all'
    elif args.snapshots!='all':
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]    
    units = args.units if args.units else 'code'
    los = args.los if args.los else '0'
    if los.lower()=='x': los = 0
    elif los.lower()=='y': los = 1
    elif los.lower()=='z': los = 2
    elif los in {'0','1','2'}:
        los = int(los)
    else: raise ValueError("Line of sight (los) must be 0, 1, 2, 'x', 'y', 'z', 'X', 'Y', 'Z'")
    c = args.c if args.c else 1
    L_0 = args.L_0 if args.L_0 else 'code'
    t_0 = args.t_0 if args.t_0 else 'Myr'

    # Load class and global variables.
    pd = PlutoData(Path)
    
    slc = args.slc if args.slc else pd.shape[los]//2

    x, y = pd.get_x('x',L_0)[::c], pd.get_x('y',L_0)[::c]
    labels = pd.get_times(snapshots,t_0)
    labels = [str(round(label,1))+t_0 for label in labels]

    savePath = os.path.join(Path,'visualise')
    if not os.path.exists(savePath):
        os.mkdir(savePath)

    plot_raw()
    plot_derived()

    