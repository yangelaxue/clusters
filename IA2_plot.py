"""
Script that plots central slices of ENZO Itasca IA2 raw fields or derived fields.
TODO: Does not take into account comoving coodinates.

Author: Angela Xue
Date: October 2026
"""

def calc_vminsvmaxs(ia2,varNames,snapshots):

    vmins, vmaxs = {}, {}

    for varName in varNames:
        if varName in {"x-velocity","y-velocity","z-velocity","Bx","By","Bz"}:
            vals = ia2.gen_vals(varName=varName,snapshots=snapshots,units=units,c=c)
            val_max = 0.
            for val in vals:
                val_max = max(val_max,np.abs(val).max())
            val_min = -val_max
        else:
            vals = ia2.gen_vals(varName=varName,snapshots=snapshots,units=units,c=c)
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
    parser.add_argument("cluster", help="Which ENZO Itasca IA2 cluster to save the properties for.", type=str)
    parser.add_argument("-DEBUG", help="Switch to debug mode.",type=bool)
    parser.add_argument("-c",type=int)
    args = parser.parse_args()
    cluster = args.cluster
    DEBUG = args.DEBUG

    c = 1
    Delta = 200

    import numpy as np
    import matplotlib.pyplot as plt
    import os

    from utils.ia2 import IA2, IA2Data
    from utils.visualise import Style, plot_vars, get_slice
    Style.set()

    ia2 = IA2Data(cluster)
    c = 1
    snapshots = IA2.snapshots[ia2.cluster]
    if DEBUG:
        c = args.c if args.c else 40
        snapshots = [10,11]
    if not os.path.exists(os.path.join(ia2.Path,"images")):
        os.mkdir(os.path.join(ia2.Path,"images"))

    # snapshots = ia2.snapshots
    units = "cgs"
    redshifts = [IA2.get_redshift(s) for s in snapshots]
    centres = [ia2.get_centre(s,c=c) for s in snapshots]
    loss = [0,1,2]

    varNames = ["Bx","By","Bz","x-velocity","y-velocity","z-velocity","Density","Dark_Matter_Density","Temperature"]
    vmins, vmaxs = calc_vminsvmaxs(ia2,varNames,snapshots)
    vmins["Density"] = None
    vmins["Dark_Matter_Density"] = None
    vmins["Temperature"] = None
    for los in loss:
        for i,s in enumerate(snapshots):
            vals = ia2.gen_vals(varNames=varNames,snapshot=s,units=units,c=c)
            vals = [get_slice(val,slc=centres[i][los],los=los) for val in vals]
            fig = plot_vars(vals,varNames,vmins=vmins,vmaxs=vmaxs)
            fig.suptitle(f"z={redshifts[i]}")
            plt.savefig(os.path.join(ia2.Path,"images",f"fields{los}_{s}.png"),dpi=Style.figkwargs["dpi"])
            plt.close()

    derived_varNames = ["ek","vorticity2","B2","curlB2"]
    vmins, vmaxs = calc_vminsvmaxs(ia2,derived_varNames,snapshots)
    for los in loss:
        for i,s in enumerate(snapshots):
            vals = ia2.gen_vals(varNames=derived_varNames,snapshot=s,units=units,c=c)
            vals = [get_slice(val,slc=centres[i][los],los=los) for val in vals]
            fig = plot_vars(vals,derived_varNames,vmins=vmins,vmaxs=vmaxs)
            fig.suptitle(f"z={redshifts[i]}")
            plt.savefig(os.path.join(ia2.Path,"images",f"derivedfields{los}_{s}.png"),dpi=Style.figkwargs["dpi"])
            plt.close()