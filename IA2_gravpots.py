"""
Script that calculates the gravitational potential of each ENZO Itasca IA2 cluster at specified snapshots.
This is done at the full resolution of the IA2 clusters.
The gravitational potential is given in CGS units.

Author: Angela Xue
Date: October 2026
"""

import numpy as np

def extend_domain(shape_ext,val):
    """
    Extend domain to desired shape by extending opposite sides equally.
    """
    
    pads = tuple(sh_bg-sh for sh,sh_bg in zip(val.shape,shape_ext))
    slices = tuple(slice(pad//2,-pad//2) for pad in pads)
    slices = []
    for pad,sh in zip(pads,shape_ext):
        if pad==0:
            slices.append(slice(0,sh))
        else:
            slices.append(slice(pad//2,-pad//2))
    slices = tuple(slices)
    ret = np.zeros(shape_ext)
    ret[slices] += val

    return ret

if __name__=="__main__":

    import os

    from utils.enzo import IA2, IA2Data, get_redshift
    from utils.gravity import get_gravpot
    from utils.units import CGS

    #%% Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("cluster", help="Which IA2 cluster to process.",type=str)
    parser.add_argument("snapshots", help="Which snapshots to animate.")
    args = parser.parse_args()
    cluster = args.cluster
    snapshots = [int(s) for s in args.snapshots.removeprefix('[').removesuffix(']').split(',')]

    ia2 = IA2Data(cluster)

    rho_gen = ia2.gen_vals("Density",snapshots,units="cgs")
    dm_gen = ia2.gen_vals("Dark_Matter_Density",snapshots,units="cgs")

    for s in snapshots:
        z = get_redshift(s)
        dL = (IA2.dL*1e3*CGS.pc) / (1+z)
        dens = next(rho_gen) + next(dm_gen)
        pot = get_gravpot(extend_domain(tuple(sh*2 for sh in dens.shape),dens),(dL,)*3,G=CGS.G)

        np.savetxt(os.path.join(ia2.Path,f"pot{s}.txt"),pot.flatten())