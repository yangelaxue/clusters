"""
Script that plots field values and derived field values from PLUTO outputs.

Author: Angela Xue
Date: October 2026
"""

def calc_vminsvmaxs(pd,varNames,snapshots,units="code"):

    vmins, vmaxs = {}, {}

    for varName in varNames:
        if varName in {"vx1","vx2","vx3","Bx1","Bx2","Bx3"}:
            vals = pd.gen_vals(varName=varName,snapshots=snapshots,units=units)
            val_max = 0.
            for val in vals:
                val_max = max(val_max,np.abs(val).max())
            val_min = -val_max
        else:
            vals = pd.gen_vals(varName=varName,snapshots=snapshots,units=units)
            val_max = 0.
            val_min = 0.
            for val in vals:
                val_max = max(val_max,val.max())
                val_min = min(val_min,val.min())
        vmins[varName] = val_min
        vmaxs[varName] = val_max
    return vmins, vmaxs

if __name__=="__main__":

    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Where the PLUTO outputs are stored.", type=str)
    args = parser.parse_args()
    Path = args.Path

    import numpy as np
    import matplotlib.pyplot as plt
    import os

    from utils.pluto import PlutoData
    from utils.visualise import Style, plot_vars
    Style.set()

    if not os.path.exists(os.path.join(Path,"images")):
        os.mkdir(os.path.join(Path,"images"))

    pd = PlutoData(Path)
    snapshots = pd.get_snapshots().tolist()
    times = pd.get_times(snapshots)
    
    varNames = ["Bx1","Bx2","Bx3","vx1","vx2","vx3","rho","prs","pot"]
    vmins, vmaxs = calc_vminsvmaxs(pd,varNames,snapshots,units="code")
    vmins["rho"] = None
    vmins["prs"] = None
    for s in snapshots:
        vals = pd.gen_vals(varNames=varNames,snapshot=s,units="code")
        fig = plot_vars(vals,varNames,vmins=vmins,vmaxs=vmaxs)
        fig.suptitle(f"t={times[s]}")
        plt.savefig(os.path.join(Path,"images",f"fields_{s}.png"),dpi=Style.figkwargs["dpi"])
        plt.close()

    derived_varNames = ["ek","vorticity2","B2","curlB2"]
    vmins, vmaxs = calc_vminsvmaxs(pd,derived_varNames,snapshots,units="code")
    for s in snapshots:
        vals = pd.gen_vals(varNames=derived_varNames,snapshot=s,units="code")
        fig = plot_vars(vals,derived_varNames,vmins=vmins,vmaxs=vmaxs)
        fig.suptitle(f"t={times[s]}")
        plt.savefig(os.path.join(Path,"images",f"derivedfields_{s}.png"),dpi=Style.figkwargs["dpi"])
        plt.close()
