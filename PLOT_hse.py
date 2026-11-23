"""
Plot the HSE status of clusters.

Author: Angela Xue
Date: October 2026
"""

import numpy as np
import os

def load_struct(fName,shape):
    import struct
    n_points = np.prod(shape)
    with open(fName,'rb') as f_o:
        val = np.array(struct.unpack('<'+'d'*n_points,f_o.read())).reshape(shape).T
    return val

def read_information(Path):
    """ Reads the topluto.info file to determine parameters used to convert raw data to PLUTO ICs. """
    global cluster, snapshots
    global rho_0, L_0, v_0
    global t_0, p_0, B_0
    global L, shape
    with open(os.path.join(Path,'topluto.info'),'r') as info:
        for line in info.readlines():
            line = line.split('#')[0].replace(" ","").replace("\n","").replace("\t","")
            if line.startswith('cluster'):
                cluster = line.split('=')[-1]
                continue
            if line.startswith('snapshots'):
                snapshots = [int(s) for s in line.split("=")[-1].removeprefix('[').removesuffix(']').split(',')]
                continue
    with open(os.path.join(Path,'info.txt'),'r') as info:
        for line in info.readlines():
            line = line.split('#')[0].replace("\n","").replace("\t","")
            if line.startswith('rho_0'):
                rho_0 = float(line.split('=')[-1])
                continue
            if line.startswith('L_0'):
                L_0 = float(line.split('=')[-1])
                continue
            if line.startswith('v_0'):
                v_0 = float(line.split('=')[-1])
                continue
            if line.startswith('shape'):
                line = line.split('=')[-1].replace(" ","")
                shape = tuple(int(sh) for sh in line.removeprefix('(').removesuffix(')').split(','))
                continue
            if line.startswith('L ='):
                line = line.split(',')[0].split("=")[-1]
                L = float(line.removeprefix("[").removesuffix("]").split(" ")[2])
                continue
    t_0 = L_0/v_0
    p_0 = rho_0*v_0**2
    B_0 = v_0 * (4*np.pi*rho_0)**.5

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Where the PLUTO run is stored.", type=str)
    parser.add_argument("-normalise", help="Divide the HSE equation by density.", type=bool)
    parser.add_argument("-input", help="Whether to plot the PLUTO initial conditions.", type=bool)
    parser.add_argument("-output", help="Whether to plot the PLUTO outputs.", type=bool)
    args = parser.parse_args()
    Path = args.Path
    normalise = args.normalise
    input = args.input
    output = args.output

    normalise = normalise if normalise else False

    if (not input) and (not output):
        print("Neither the PLUTO initial conditions not outputs are to be examined.")
        exit()

    if not os.path.exists(os.path.join(Path,"images")):
        os.mkdir(os.path.join(Path,"images"))

    from utils.units import CGS

    
    

    

    # import matplotlib.pyplot as plt
    # import os

    from utils.gradient import gradient_radial
    from utils.gravity import calc_gravforce_from_pot

    # from utils.units import CGS
    from utils.grid import get_xyz

    from utils.visualise import Style, plot_var, get_slice
    Style.set()


    if input:
        read_information(Path)
        dL = L/shape[0]
        centre = tuple(sh//2 for sh in shape)
    
        rho = load_struct(fName=os.path.join(Path,"rho0.dbl"),shape=shape)
        prs = load_struct(fName=os.path.join(Path,"prs0.dbl"),shape=shape)
        pot = load_struct(fName=os.path.join(Path,"pot00.dbl"),shape=shape)

        LHS = gradient_radial(prs,centre,dxdydz=dL)
        grav = np.sum([gravx**2 for gravx in calc_gravforce_from_pot(pot,dL)],axis=0)**.5
        RHS = -grav*rho

        if normalise:
            title = r"$(\nabla P-(-g\rho))/\rho$"
            val = get_slice((LHS-RHS)/rho,slc=64,los="y")
        else:
            title = r"$\nabla P-(-g\rho)$"
            val = get_slice((LHS-RHS),slc=64,los="y")
        vmax = np.max(val[np.abs(val)<np.inf])
        vmin = -vmax
        x, y = get_xyz(shape[:2])
        x = x*dL
        y = y*dL

        fig = plot_var(val,"",cmap="bwr",vmin=vmin,vmax=vmax,title=title,x=x,y=y)
        fig.savefig(os.path.join(Path,"images","HSE_input.png"))
    if output:
        from utils.pluto import PlutoData
        pd = PlutoData(Path)
        c = 1
        units = "code"
        dL = pd.get_dxdydz()[0]
        centre = tuple(sh//2 for sh in pd.shape)
        dxdydz = pd.get_dxdydz(units=units)[0]

        snapshots = pd.get_snapshots()
        rhos = pd.gen_vals(varName="rho",snapshots=snapshots,units=units)
        prss = pd.gen_vals(varName="prs",snapshots=snapshots,units=units)
        pots = pd.gen_vals(varName="pot",snapshots=snapshots,units=units)
        vals = []

        for rho,prs,pot in zip(rhos,prss,pots):
            LHS = gradient_radial(prs,centre,dxdydz=dxdydz)
            grav = np.sum([gravx**2 for gravx in calc_gravforce_from_pot(pot,dxdydz)],axis=0)**.5
            RHS = -grav*rho
            if normalise:
                title = r"$(\nabla P-(-g\rho))/\rho$"
                val = get_slice((LHS-RHS)/rho,slc=64,los="y")
            else:
                title = r"$\nabla P-(-g\rho)$"
                val = get_slice((LHS-RHS),slc=64,los="y")
            vals.append(val)

        vmax = np.abs(vals).max()
        vmin = -vmax
        x, y = get_xyz(pd.shape[:2])
        x = x*dL
        y = y*dL
        for s,val in zip(snapshots,vals):
            fig = plot_var(val,"",cmap="bwr",vmin=vmin,vmax=vmax,title=title,x=x,y=y)
            fig.savefig(os.path.join(Path,"images",f"HSE_output_{s}.png"))
