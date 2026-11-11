"""
Script that saves central slices of ENZO Itasca IA2 raw fields or derived fields.

Author: Angela Xue
Date: October 2026
"""

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
    import os, pickle

    from utils.ia2 import IA2, IA2Data
    from utils.visualise import Style, get_slice
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
    for los in loss:
        for i,s in enumerate(snapshots):
            vals = ia2.gen_vals(varNames=varNames,snapshot=s,units=units,c=c)
            vals = [get_slice(val,slc=centres[i][los],los=los) for val in vals]
            for varName,val in zip(varNames,vals):
                pickle.dump(val,open(os.path.join(ia2.Path,"images",f"{varName}_slice{los}_{s}.p"),"wb"))

    derived_varNames = ["ek","vorticity2","B2","curlB2"]
    for los in loss:
        for i,s in enumerate(snapshots):
            vals = ia2.gen_vals(varNames=derived_varNames,snapshot=s,units=units,c=c)
            vals = [get_slice(val,slc=centres[i][los],los=los) for val in vals]
            for varName,val in zip(derived_varNames,vals):
                pickle.dump(val,open(os.path.join(ia2.Path,"images",f"{varName}_slice{los}_{s}.p"),"wb"))