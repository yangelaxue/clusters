"""
Script to calculate the rms velocity and/or Mach number of a whole field, then plot.

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
    rhos = pd.gen_vals('rho',snapshots)
    prss = pd.gen_vals('prs',snapshots)
    gamma = pd.gamma
    Mach_rms = np.array([get_mach_rms(vx1,vx2,vx3,subtract_mean=True,normalise=calcMach,rho=rho,prs=prs,gamma=gamma)[0] for vx1,vx2,vx3,rho,prs in zip(vx1s,vx2s,vx3s,rhos,prss)])
    return Mach_rms

if __name__=="__main__":

    import numpy as np
    from matplotlib import pyplot as plt
    import os, argparse

    from utils.pluto import PlutoData
    from utils.visualize import DefaultStyle, set_style
    from utils.fluids import get_mach_rms, get_val_rms

    set_style()
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("-snapshots", help="Which snapshots to animate.",required=False)
    parser.add_argument("-subtract_mean", help="Whether to subtract mean before calculating value RMS. Defaults to False.",required=False)

    args = parser.parse_args()

    # Validate inputs.
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path

    pd = PlutoData(Path)

    if not args.snapshots:
        snapshots = [i for i in range(len(pd.times))]
    else:
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]
    subtract_mean = args.subtract_mean if args.subtract_mean else False

    savePath = os.path.join(Path,'visualise')
    if not os.path.exists(savePath):
        os.mkdir(savePath)

    times = np.array([pd.times[s] for s in snapshots])

    # Calculate Mach_rms
    Mach_rms = get_vrms(calcMach=True)

    # Calculate rho_rms
    rho = pd.gen_vals('rho',snapshots=snapshots)
    rho_rms = np.array([get_val_rms(_rho,subtract_mean,normalise=True) for _rho in rho])

    # Calculate prs_rms
    prs = pd.gen_vals('prs',snapshots=snapshots)
    prs_rms = np.array([get_val_rms(_prs,subtract_mean,normalise=True) for _prs in prs])

    # Save to txt file.
    arr = np.concatenate([times[:,np.newaxis],Mach_rms[:,np.newaxis],rho_rms[:,np.newaxis],prs_rms[:,np.newaxis]],axis=1)
    np.savetxt('rms.txt',arr,header='times Mach_rms rho_rms prs_rms')

    