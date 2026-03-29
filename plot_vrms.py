"""
Script to calculate the rms velocity and/or Mach number of a whole field, then save and plot.

Author: Angela Xue
Date: March 2026
"""

def get_vrms(calcMach=False):
    """
    Returns an ordered list of root-mean-square velocity or turbulent Mach number
    """

    vx1s = pd.gen_vals('vx1',snapshots)
    vx2s = pd.gen_vals('vx2',snapshots)
    vx3s = pd.gen_vals('vx3',snapshots)
    if calcMach:
        rhos = pd.gen_vals('rho',snapshots)
        prss = pd.gen_vals('prs',snapshots)
        gamma = pd.gamma
        ret = np.array([get_mach_rms(vx1,vx2,vx3,subtract_mean=True,normalise=calcMach,rho=rho,prs=prs,gamma=gamma) for vx1,vx2,vx3,rho,prs in zip(vx1s,vx2s,vx3s,rhos,prss)])
        Machrms, c_s = tuple(ret.transpose().tolist())
        return np.array(Machrms), np.array(c_s)
    else:
        return np.array([get_mach_rms(vx1,vx2,vx3,subtract_mean=True,normalise=calcMach) for vx1,vx2,vx3 in zip(vx1s,vx2s,vx3s)])

if __name__=="__main__":

    import numpy as np
    from matplotlib import pyplot as plt
    import os, argparse

    from utils.pluto import PlutoData
    from utils.visualize import DefaultStyle, set_style
    from utils.fluids import get_mach_rms

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
    Mach_rms, c_s = get_vrms(True)
    v_rms = Mach_rms * c_s

    arr = np.array([pd.get_times(snapshots), v_rms, c_s, Mach_rms]).T

    # Save to text file
    np.savetxt(os.path.join(pd.Path,'vrms.txt'),arr,header='time v_rms c_s Mach_rms')

    # Plot root-mean-square velocity
    fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'],tight_layout=True)

    ax.plot(times,v_rms,marker='x')
    ax.plot(times,c_s,label=r'$c_s$',marker='^')

    ax.grid()
    ax.set_xlabel(f'Time ({t_0})')
    ax.set_ylabel('Speed (km/s)')
    ax.set_title('Root-mean-square Velocity')
    ax.legend(frameon=False)

    plt.savefig(os.path.join(savePath,f'vrms_test.png'),dpi=DefaultStyle.figkwargs['dpi'])

    # Plot turbulent Mach number
    fig, ax = plt.subplots(figsize=DefaultStyle.figkwargs['figsize'],tight_layout=True)

    ax.plot(times,Mach_rms,marker='x')

    ax.grid()
    ax.set_xlabel(f'Time ({t_0})')
    ax.set_ylabel('Mach Number')
    ax.set_title('Root-mean-square Mach Number')

    plt.savefig(os.path.join(savePath,f'machrms_test.png'),dpi=DefaultStyle.figkwargs['dpi'])

