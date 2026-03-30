"""
Inspired by Mohapatra & Sharma (2019), create plots power spectra.

Author: Angela Xue
Date: Oct 2025
"""

def plot(P,varName):

    fig, ax = plt.subplots(**DefaultStyle.figkwargs)

    if varName=='rho':
        label = r"\rho"
    elif varName=='prs':
        label = r"P"
    elif varName=='v':
        label = r"v"
    elif varName=='xraySB':
        label = r"{\rm SB}_{\rm X}"
    elif varName=='SZSB':
        label = r"{\rm SB}_{\rm SZ}"
    elif varName=='eta':
        pass
    else:
        raise NotImplementedError(f"Variable name {varName} not implemented.")

    if subtract_mean:
        label = r"\delta "+label
    if normalise:
        label = label + fr"/\langle {label}\rangle"
    
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(fr"$P_{label}(k)k^2$")
    if varName=='eta':
        ax.set_ylabel(fr"$\eta(k)$")

    for i,_P in enumerate(P):
        alpha = ((i+1)/len(snapshots))**1
        ax.plot(k,_P,c=colors(alpha))

    plt.savefig(os.path.join(savePath,f'powerspectrum_{varName}.png'),dpi=DefaultStyle.figkwargs['dpi'])
    plt.close()

def save(P,varName):

    arr = np.concatenate([times[:,np.newaxis],P],axis=1)
    np.savetxt(os.path.join(pd.Path,f'powerspectrum_{varName}.txt'),arr,header='k P_k')

if __name__=="__main__":

    import numpy as np
    from matplotlib import pyplot as plt
    import os, argparse

    from utils.fluids import get_c_s
    from utils.pluto import PlutoData
    from utils.visualize import DefaultStyle, set_style

    from Kea.statistics.spectra.spectra_base import calculate_integrated_spectrum

    set_style()
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("-snapshots", help="Which snapshots to plot.",required=False)
    parser.add_argument("-L_0", help="Length unit to plot.", type=int,required=False)
    parser.add_argument("-subtract_mean", help="Whether to subtract mean val.", type=int,required=False)
    parser.add_argument("-normalise", help="Whether to normalise spectrum.", type=int,required=False)

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
    L_0 = args.L_0 if args.L_0 else 'kpc'
    subtract_mean = int(args.subtract_mean) if args.subtract_mean else 0
    normalise = int(args.normalise) if args.normalise else 1
    
    savePath = os.path.join(Path,'visualise')
    if not os.path.exists(savePath):
        os.mkdir(savePath)

    times = np.array([pd.times[s] for s in snapshots])
    lenn = np.diff(pd.bbox,axis=1).flatten()
    k = calculate_integrated_spectrum([np.ones(pd.shape)],lenn=lenn)[0]
    colors = plt.cm.YlGnBu

    # Plot and save density power spectrum
    rho = pd.gen_vals('rho',snapshots=snapshots)
    if normalise:
        P_rho = np.array([calculate_integrated_spectrum([_val],lenn=lenn)[1]/_val.mean()**2 for _val in rho])
    else:
        P_rho = np.array([calculate_integrated_spectrum([_val],lenn=lenn)[1] for _val in rho])
    plot(P_rho,'rho')
    save(P_rho,'rho')

    # Plot and save pressure power spectrum
    if 'prs' in pd.varNames:
        prs = pd.gen_vals('prs',snapshots=snapshots)
        if normalise:
            P = np.array([calculate_integrated_spectrum([_val],lenn=lenn)[1]/_val.mean()**2 for _val in prs])
        else:
            P = np.array([calculate_integrated_spectrum([_val],lenn=lenn)[1] for _val in prs])
        plot(P,'prs')
        save(P,'prs')

    # Plot and save velocity power spectrum
    vx1 = pd.gen_vals('vx1',snapshots=snapshots)
    vx2 = pd.gen_vals('vx2',snapshots=snapshots)
    if pd.ndim==2:
        v = ((_vx1**2 + _vx2**2)**.5 for _vx1,_vx2 in zip(vx1,vx2))
    elif pd.ndim==3:
        vx3 = pd.gen_vals('vx3',snapshots=snapshots)
        v = ((_vx1**2 + _vx2**2 + _vx3**2)**.5 for _vx1,_vx2,_vx3 in zip(vx1,vx2,vx3))
    if normalise:
        rho = pd.gen_vals('rho',snapshots=snapshots)
        prs = pd.gen_vals('prs',snapshots=snapshots)
        c_s = [get_c_s(_rho,_prs,gamma=5/3) for _rho,_prs in zip(rho,prs)]
        P_v = np.array([calculate_integrated_spectrum([_val],lenn=lenn)[1]/_c_s**2 for _val,_c_s in zip(v,c_s)])
    else:
        P_v = np.array([calculate_integrated_spectrum([_val],lenn=lenn)[1] for _val in v])
    plot(P_v,'v')
    save(P_v,'v')

    plot(P_rho/P_v,'eta')


    
