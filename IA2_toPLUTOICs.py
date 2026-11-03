"""
Script that generates PLUTO initial conditions from ENZO Itasca IA2 clusters.

Author: Angela Xue
Date: October 2026
"""

import numpy as np
from matplotlib import pyplot as plt
import os, struct

from utils.units import CGS
from utils.grid import get_Rmesh

def read_topluto(Path):
    """ Reads the topluto.info file to determine parameters used to convert raw data to PLUTO ICs. """
    global cluster, snapshots, c, Delta, mu, apd_radius, apd_slope, L_final, shape_final
    cluster, snapshots, c, Delta, mu, apd_radius, apd_slope, L_final, shape_final = \
        "E14", "10,11", 40, 200, 0.5, 0.5, 16., 3.0, (64,64,64)
    with open(os.path.join(Path,"topluto.info"),"r") as info:
        for line in info.readlines():
            line = line.split("#")[0].replace(" ","").replace("\n","").replace("\t","")
            if line.startswith("cluster"):
                cluster = line.split("=")[-1]
                continue
            if line.startswith("snapshots"):
                snapshots = [int(s) for s in line.split("=")[-1].removeprefix("[").removesuffix("]").split(",")]
                continue
            if line.startswith("c"):
                c = int(line.split("=")[-1])
                continue
            if line.startswith("Delta"):
                Delta = int(line.split("=")[-1])
                continue
            if line.startswith("mu"):
                mu = float(line.split("=")[-1]) # Mean molecular mass
                continue
            if line.startswith("L_final"):
                L_final = float(line.split("=")[-1])
                continue
            if line.startswith("shape_final"):
                line = line.split("=")[-1]
                try:
                    shape_final = (int(line),)*3
                except:
                    shape_final = tuple(int(sh) for sh in line.removeprefix("(").removesuffix(")").split(","))
                continue
            if line.startswith("apd_radius"):
                apd_radius = float(line.split("=")[-1])
                continue
            if line.startswith("apd_slope"):
                apd_slope = float(line.split("=")[-1])
                continue

def get_shapes_ext(ia2,snapshots,c):
    redshifts = tuple(IA2.get_redshift(s) for s in snapshots)
    centres = tuple(ia2.get_centre(s,c) for s in snapshots)
    L_max = []
    for z,cen in zip(redshifts,centres):
        L = ia2.dim*IA2.dL / (1+z) / L_0
        cen = tuple(_cen*c*IA2.dL / (1+z) / L_0 for _cen in cen)
        L_max.append(max(tuple(max(_cen,L-_cen) for _cen in cen)) * 2)
    L_max = max(L_max) # Cluster units.
    shapes_ext = tuple((int(L_max*(1+z)/(IA2.dL/L_0))//c,)*3 for z in redshifts)
    return shapes_ext

def apodize(val,rad,sl,dL):
    """ Apodize a cluster about the grid centre. """
    R = get_Rmesh(val.shape)*dL
    ret = val * 1/(1+np.exp(sl*(R-rad)))
    return ret

def fade_outskirts(val,dL):
    """ Fade the outskirts of a cluster to a beta profile. """
    from scipy.optimize import curve_fit
    from utils.clusters import beta_profile
    centre = tuple(idx[0] for idx in np.where(val==val.max()))
    R = get_Rmesh(shapes_ext[0],centre)*dL
    beta_args = curve_fit(beta_profile,R[val>0].flatten(),val[val>0].flatten(),p0=(val.max(),3,1.5))[0]
    profile = beta_profile(R,1.,*beta_args[1:])
    profile[profile>1e-2] = 1e-2
    profile = profile * beta_args[0]
    ret = apodize(val,1.0,16.,dL)
    ret[ret<=profile] = profile[ret<=profile]
    return ret

def interpolate(val,L_in,L_out,shape_out):
    """ Interpolate field to desired shape. """
    from scipy.interpolate import RegularGridInterpolator
    xyz_in = tuple(np.linspace(-L_in/2,L_in/2,sh) for sh in val.shape)
    xyz_out = tuple(np.linspace(-L_out/2,L_out/2,sh) for sh in shape_out)
    points = (np.linspace(_x[0],_x[-1],sh) for _x,sh in zip(xyz_out,shape_out))
    points = np.vstack([X.ravel() for X in np.meshgrid(*points,indexing="ij")]).T
    ret = RegularGridInterpolator(xyz_in,val)(points).reshape(shape_out)
    return ret

def plot(los):
    """ Plot the final apodised fields. """
    from matplotlib.colors import LogNorm
    from utils.visualise import Style, get_slice
    Style.set()

    plotNames = ["vx1","vx2","vx3","Bx1","Bx2","Bx3","rho","prs"]
    for i in range(len(snapshots)):
        plotNames.append(f"pot{i}")
    def load_vals():
        for varName in plotNames:
            with open(os.path.join(savePath, f"{varName}0.dbl"),"rb") as f_o:
                val = np.array(struct.unpack("<"+"d"*n_points,f_o.read())).reshape(shape_final).T
            yield val
    vals = load_vals()

    nrows = int(np.ceil(len(plotNames)/3))
    fig = plt.figure(figsize=(10,nrows*3),tight_layout=True)
    for i,(varName,val) in enumerate(zip(plotNames,vals)):
        varkwargs = Style.get_varkwargs(varName)
        ax = fig.add_subplot(nrows,3,i+1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(varkwargs["title"])
        val = get_slice(val,slc=val.shape[0]//2,los=los,c=1)
        norm = None
        vmin, vmax = None, None
        if varName in {"rho","dm","prs"}: #TODO
            norm = LogNorm()
        if varName.startswith("vx") or varName.startswith("Bx"):
            vmax = np.abs(val).max()
            vmin = -vmax

        cmap = varkwargs["cmap"]
        im = ax.pcolormesh(val,cmap=cmap,vmin=vmin,vmax=vmax,norm=norm)
        cbar = Style.get_cbar(im,fig,ax,size=.1)

    if not os.path.exists(os.path.join(savePath,"images")):
        os.mkdir(os.path.join(savePath,"images"))

    plt.savefig(os.path.join(savePath,"images",f"cluster_{los}.png"),dpi=300)
    plt.close()

if __name__=="__main__":
    
    import os, struct

    from utils.ia2 import IA2, IA2Data
    from utils.clusters import calc_avgvelocity
    from utils.gravity import calc_gravpot
    from utils.grid import extend_domain

    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("savePath", help="Where to save ICs.",type=str)
    args = parser.parse_args()
    savePath = args.savePath

    # Load interpolation arguments from topluto.info
    read_topluto(savePath)
    # Load cluster.
    ia2 = IA2Data(cluster)
    # Load cluster units.
    units = "Delta"+str(Delta)
    rho_0, L_0, v_0, t_0, p_0, B_0 = ia2.get_units(Delta)
    K_0 = CGS.mp * v_0**2 / CGS.kB
    G_0 = CGS.G / (1/rho_0/t_0**2)

    # Load time varying cluster properties
    centres = tuple(ia2.get_centre(s,c=c) for s in snapshots)
    dLs = tuple(IA2.dL * c / (1+IA2.get_redshift(s)) / L_0 for s in snapshots)

    # Generate raw field values.
    scalarNames = ["Density","Pressure"]
    vxNames = ["x-velocity","y-velocity","z-velocity"]
    BxNames = ["Bx","By","Bz"]
    scalarvals = ia2.gen_vals(varNames=scalarNames,snapshot=snapshots[0],c=c,units=units,mu=mu)
    vxvals = ia2.gen_vals(varNames=vxNames,snapshot=snapshots[0],c=c,units=units)
    Bxvals = ia2.gen_vals(varNames=BxNames,snapshot=snapshots[0],c=c,units=units)
    dens = (
        rho+dm for rho,dm in zip(
            ia2.gen_vals("Density",snapshots=snapshots,c=c,units=units),
            ia2.gen_vals("Dark_Matter_Density",snapshots=snapshots,c=c,units=units)
            )
        )
    rho = ia2.get_val("Density",snapshot=snapshots[0],c=c,units=units)

    # Extend volume to the physical dimensions.
    shapes_ext = get_shapes_ext(ia2,snapshots,c=c)
    scalarvals = (extend_domain(val,shapes_ext[0],centres[0]) for val in scalarvals)
    vxvals = (extend_domain(val,shapes_ext[0],centres[0]) for val in vxvals)
    Bxvals = (extend_domain(val,shapes_ext[0],centres[0]) for val in Bxvals)
    dens = (extend_domain(val,shapes_ext[i],centres[i]) for i,val in enumerate(dens))
    rho = extend_domain(rho,shapes_ext[0],centres[0])
    
    # Remove average velocity.
    vxvals = (vx-calc_avgvelocity(rho,vx) for vx in vxvals)
    # Calculate gravitational potentials.
    pot = (calc_gravpot(dens,(dL,)*3,G=G_0) for dens,dL in zip(dens,dLs))
    pot = (pot-pot.max() for pot in pot)

    # Apodize all fields
    scalarvals = (fade_outskirts(val,dLs[0]) for val in scalarvals)
    if apd_radius>0:
        vxvals = (apodize(val,apd_radius,apd_slope,dLs[0]) for val in vxvals)
        Bxvals = (apodize(val,apd_radius,apd_slope,dLs[0]) for val in Bxvals)

    # Interpolate to final dimensions.
    Ls = tuple(shapes_ext[i][0]*dLs[i] for i in range(len(snapshots)))
    L_final = min(L_final,*Ls)
    scalarvals = (interpolate(val,Ls[0],L_final,shape_final) for val in scalarvals)
    Bxvals = (interpolate(val,Ls[0],L_final,shape_final) for val in Bxvals)
    vxvals = (interpolate(val,Ls[0],L_final,shape_final) for val in vxvals)
    pot = (interpolate(val,L,L_final,shape_final) for val,L in zip(pot,Ls))

    # Save final fields!
    # Find critical density.
    # H2_a = lambda a : IA2.H_0**2 * (IA2.Omega_r0/a**4 + (IA2.Omega_bm0+IA2.Omega_dm0)/a**3 + IA2.Omega_Lambda0 + IA2.Omega_k0/a**2)
    # H2_z = lambda z : H2_a(1/(1+z))
    # rho_cr_z = lambda z: 3*H2_z(z)/(8*np.pi*CGS.G)
    # temp = ia2.get_val("Temperature",snapshot=snapshots[0],c=c)
    # rho_cr = rho_cr_z(IA2.get_redshift(snapshots[0]))/rho_0 / 2
    # prs_cr = rho_cr*rho_0/CGS.me*CGS.kB*(temp[temp>0].max()*1e-4)/p_0
    # rho_cr, prs_cr = 0., 0.

    n_points = np.prod(shape_final)
    for varName,val in zip(scalarNames,scalarvals):
        print("Processing variable:",varName)
        if varName=="Density":
            varName = "rho"
            rho_min = val.min()
        if varName=="Pressure":
            varName = "prs"
            prs_min = val.min()
        with open(os.path.join(savePath, f"{varName}0.dbl"),"wb") as f_o:
            f_o.write(struct.pack("<"+"d"*n_points,*(val.T.flatten())))
    for varName,val in zip(vxNames,vxvals):
        print("Processing variable:",varName)
        if varName.startswith("x"):
            varName = "vx1"
        elif varName.startswith("y"):
            varName = "vx2"
        elif varName.startswith("z"):
            varName = "vx3"
        with open(os.path.join(savePath, f"{varName}0.dbl"),"wb") as f_o:
            f_o.write(struct.pack("<"+"d"*n_points,*(val.T.flatten())))
    for varName,val in zip(BxNames,Bxvals):
        print("Processing variable:",varName)
        if varName.endswith("x"):
            varName = "Bx1"
        elif varName.endswith("y"):
            varName = "Bx2"
        elif varName.endswith("z"):
            varName = "Bx3"
        with open(os.path.join(savePath, f"{varName}0.dbl"),"wb") as f_o:
            f_o.write(struct.pack("<"+"d"*n_points,*(val.T.flatten())))
    # Save potentials velocities
    for i,val in enumerate(pot):
        print("Processing variable:","pot"+str(i))
        with open(os.path.join(savePath, f"pot{i}0.dbl"),"wb") as f_o:
            f_o.write(struct.pack("<"+"d"*n_points,*(val.T.flatten())))
    # Save grid.
    DIM = 3
    L = (L_final,)*3
    with open(os.path.join(savePath, "grid0.out"),"w") as f_o:
        f_o.write("# GEOMETRY:   CARTESIAN\n")
        for d in range(DIM):
            f_o.write(f"{shape_final[d]}\n")
            for i in range(shape_final[d]):
                xL = -0.5 + (i-0.5)/(shape_final[d] - 1.)
                xR = -0.5 + (i+0.5)/(shape_final[d] - 1.)
                xL *= L[d]
                xR *= L[d]
                f_o.write("{:d}   {:12.6e}  {:12.6e}\n".format(i+1, xL, xR))

    def get_times():
        from astropy.cosmology import FlatLambdaCDM
        import astropy.units as u
        cosmo = FlatLambdaCDM(H0=IA2.h_0*100*u.km/u.s/u.Mpc, Om0=IA2.Omega_dm0+IA2.Omega_bm0)

        times = []
        for s in snapshots:
            z = IA2.get_redshift(s)
            times.append(cosmo.age(z).value*(1e9*CGS.yr))
        times = np.array(times)
        return times-times.min()
    times = get_times()
    # prs_min = rho_cr/CGS.mp * CGS.kB * 2.75
    with open(os.path.join(savePath,"info.txt"),"w") as info:
        info.write(f"L = {np.array(L)}, L/2 = {np.array(L)/2}\n")
        info.write(f"rho_0 = {rho_0}\n")
        info.write(f"L_0 = {L_0}\n")
        info.write(f"v_0 = {v_0}\n")
        info.write(f"shape = {shape_final}\n")
        info.write("\n")
        info.write(f"rho_min = {rho_min}\n")
        info.write(f"prs_min = {prs_min}\n")
        info.write("\n")
        info.write(f"times = {times/t_0}\n")

    plot("x")
    plot("y")
    plot("z")