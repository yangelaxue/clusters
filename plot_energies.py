"""
Script that sums the different energy components in a turbulent field.
Inspired by Tulasi Parashar's script in TurbAn.

Author: Angela Xue
Date: March 2026
"""

def get_ek():
    """ Bulk kinetic energy. """
    rho = pd.gen_vals('rho',snapshots=snapshots,units='cgs')
    vx1 = pd.gen_vals('vx1',snapshots=snapshots,units='cgs')
    vx2 = pd.gen_vals('vx2',snapshots=snapshots,units='cgs')
    if pd.ndim==2:
        v2 = (_vx1**2 + _vx2**2 for _vx1,_vx2 in zip(vx1,vx2))
    elif pd.ndim>2:
        vx3 = pd.gen_vals('vx3',snapshots=snapshots,units='cgs')
        v2 = (_vx1**2 + _vx2**2 + _vx3**2 for _vx1,_vx2,_vx3 in zip(vx1,vx2,vx3))
    return np.array([(0.5*_rho*_v2).sum()*dV for _rho,_v2 in zip(rho,v2)])

def get_et():
    """ Calculate total thermal energy. """
    prs = pd.gen_vals('prs',snapshots=snapshots,units='cgs')
    return np.array([(_prs/(pd.gamma-1)).sum()*dV for _prs in prs])

def get_eB():
    """ Calculate total magnetic energy. """
    Bx1 = pd.gen_vals('Bx1',snapshots=snapshots,units='cgs')
    Bx2 = pd.gen_vals('Bx2',snapshots=snapshots,units='cgs')
    if pd.ndim==2:
        B2 = (_Bx1**2 + _Bx2**2 for _Bx1,_Bx2 in zip(Bx1,Bx2))
    elif pd.ndim>2:
        Bx3 = pd.gen_vals('Bx3',snapshots=snapshots,units='cgs')
        B2 = (_Bx1**2 + _Bx2**2 + _Bx3**2 for _Bx1,_Bx2,_Bx3 in zip(Bx1,Bx2,Bx3))
    return np.array([(_B2/2).sum()*dV for _B2 in B2])

if __name__=="__main__":

    import numpy as np
    from matplotlib import pyplot as plt
    import os, argparse

    from utils.pluto import PlutoData
    from utils.visualize import DefaultStyle, set_style

    set_style()
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("-snapshots", help="Which snapshots to animate.",required=False)
    parser.add_argument("-t_0", help="Time unit to plot.", type=int,required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    if not args.snapshots:
        snapshots = 'all'
    else:
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]
    t_0 = args.t_0 if args.t_0 else 'Myr'

    pd = PlutoData(Path)

    savePath = os.path.join(Path,'visualise')
    if not os.path.exists(savePath):
        os.mkdir(savePath)

    times = pd.get_times(snapshots,t_0)
    dV = np.prod(np.diff(pd.bbox,axis=1).flatten()/pd.shape) * pd.L_0
    ek = get_ek()
    et = get_et()
    if pd.PHYSICS in {'MHD','RMHD'}:
        eB = get_eB()

    # Plot energies
    fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'],tight_layout=True)

    ax.plot(times,ek,label="Kinetic")
    ax.plot(times,et,label="Thermal")
    if pd.PHYSICS in {'MHD','RMHD'}:
        ax.plot(times,eB,label="Magnetic")

    ax.set_xlabel(f'Time ({t_0})')
    ax.set_ylabel('Energy (cgs)')
    ax.set_title('Energy Contributions')
    ax.grid()
    ax.legend(frameon=False)

    plt.savefig(os.path.join(savePath,f'energy.png'),dpi=DefaultStyle.figkwargs['dpi'])

    # Plot energy contributions

    etot = ek + et
    if pd.PHYSICS in {'MHD','RMHD'}:
        etot += eB
    fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'],tight_layout=True)

    ax.plot(times,ek/etot,label="Kinetic")
    ax.plot(times,et/etot,label="Thermal")
    if pd.PHYSICS in {'MHD','RMHD'}:
        ax.plot(times,eB/etot,label="Magnetic")

    ax.set_xlabel(f'Time ({t_0})')
    ax.set_title('Energy Fractions')
    ax.grid()
    ax.legend(frameon=False)

    plt.savefig(os.path.join(savePath,f'energyfrac.png'),dpi=DefaultStyle.figkwargs['dpi'])
    
    # Save energies to file.
    if pd.PHYSICS in {'MHD','RMHD'}:
        with open(os.path.join(Path,"energies.txt"),"w") as f:
            f.write("time ek ep eB\n")
            for i, time in enumerate(times):
                f.write(f"{time} {ek[i]} {et[i]} {eB[i]}")
    else:
        with open(os.path.join(Path,"energies.txt"),"w") as f:
            f.write("time ek ep\n")
            for i, time in enumerate(times):
                f.write(f"{time} {ek[i]} {et[i]}")