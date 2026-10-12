"""
Script that plots slices, projections of ENZO Itasca IA2 raw fields or derived fields.

Author: Angela Xue
Date: October 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from utils.visualise import Style, get_cbar

from utils.grid import get_xyz

def plot_slice(ia2,varName,snapshot,los):

    centre = ia2.get_centre(snapshot,c=c)
    r_Delta, M_Delta, c_s, t_s, den_centre = ia2.get_units(Delta)

    val = ia2.get_val(varName,snapshot,units="cgs",c=c)
    if varName in {"vorticity","curlB"}:
        val = np.sum([_val**2 for _val in val])
    if "Density" in varName:
        val /= den_centre
    elif varName in {"x-velocity","y-velocity","z-velocity"}:
        val /= c_s
    elif varName in {"Bx","By","Bz"}:
        pass
    xyz = get_xyz(val.shape,centre)
    if los in {0,"x","X"}:
        val = val[centre[0]]
        x, y = xyz[1], xyz[2]
    elif los in {1,"y","Y"}:
        val = val[:,centre[1]]
        x, y = xyz[2], xyz[1]
    elif los in {1,"z","Z"}:
        val = val[:,:,centre[2]]
        x, y = xyz[0], xyz[1]
    else:
        raise ValueError("Line-of-sight must be 0, x, X, 1, y, Y, 2, z, Z.")
    
    dL = IA2.dL * c / (1+IA2.get_redshift(snapshot)) / r_Delta
    x, y = x*dL, y*dL
    
    varkwargs = Style.get_varkwargs(varName)
    
    fig, ax = plt.subplots(figsize=(5,5),**Style.figkwargs)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$ $(R_{})$".format("{"+str(Delta)+"}"))
    ax.set_ylabel("$y$ $(R_{})$".format("{"+str(Delta)+"}"))
    im = ax.pcolormesh(x,y,val,cmap=varkwargs["cmap"])
    cbar = get_cbar(im,fig,ax)
    cbar.ax.set_title(r"$\rho$ $(\rho_{\rm c})$",loc='left')
    plt.show()

def plot_projection(ia2,varName,snapshot):
    pass


if __name__=="__main__":

    from utils.ia2 import IA2, IA2Data

    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("cluster", help="Which ENZO Itasca IA2 cluster to save the properties for.", type=str)
    args = parser.parse_args()
    cluster = args.cluster

    c = 1
    Delta = 200
    if 1: # DEBUG
        c = 40

    ia2 = IA2Data(cluster)
    plot_slice(ia2,"Density",15,los="x")
    # for s in ia2.snapshots:
    #     for varName in IA2.fields:
    #         plot_slice(ia2,varName,s)
    #     for varName in {"ek","vorticity","B2","curlB"}:
    #         plot_slice(ia2,varName,s)
