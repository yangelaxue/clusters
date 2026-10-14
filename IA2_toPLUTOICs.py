"""
Script that generates PLUTO initial conditions from ENZO Itasca IA2 clusters.

Author: Angela Xue
Date: October 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from utils.units import CGS

def read_topluto(Path):
    """ Reads the topluto.info file to determine parameters used to convert raw data to PLUTO ICs. """
    global c, Delta, mu, apd_radius, apd_slope, L_final, shape_final
    c, Delta, mu, apd_radius, apd_slope, L_final, shape_final = 40, 200, 0.5, 0.5, 16., 3.0, (64,64,64)
    with open(os.path.join(Path,'topluto.info'),'r') as info:
        for line in info.readlines():
            line = line.split('#')[0].replace(" ","").replace("\n","").replace("\t","")
            if line.startswith('c'):
                c = int(line.split('=')[-1])
                continue
            if line.startswith('Delta'):
                Delta = int(line.split('=')[-1])
                continue
            if line.startswith('mu'):
                mu = float(line.split('=')[-1]) # Mean molecular mass
                continue
            if line.startswith('apd_radius'):
                apd_radius = float(line.split('=')[-1])
                continue
            if line.startswith('apd_slope'):
                apd_slope = float(line.split('=')[-1])
                continue
            if line.startswith('L_final'):
                L_final = float(line.split('=')[-1])
                continue
            if line.startswith('shape_final'):
                line = line.split('=')[-1]
                try:
                    shape_final = (int(line),)*3
                except:
                    shape_final = tuple(int(sh) for sh in line.removeprefix('(').removesuffix(')').split(','))
                continue

def define_units(ia2,Delta):
    """ Define code units to use in PLUTO. These are cluster units. """
    r_Delta, M_Delta, c_s, t_s, den_centre = ia2.get_units(Delta)

    global L_0, v_0, rho_0
    global t_0, B_0, p_0
    global K_0, G_0

    L_0 = r_Delta
    v_0 = c_s
    rho_0 = den_centre
    
    t_0 = L_0/v_0 # time unit in seconds
    B_0 = v_0 * (4*np.pi*rho_0)**.5 # Magnetic field unit in Gauss (G)
    p_0 = rho_0*v_0**2 # Pressure unit in dyne/cm^2
    
    K_0 = CGS.mp * v_0**2 / CGS.kB
    G_0 = CGS.G / (1/rho_0/t_0**2)

def extend_domain(val,shape_ext,centre):
    """ Extend domain centering the cluster centre. """

    shape_init = val.shape
    center_ext = tuple(_sh//2 for _sh in shape_ext)
    slices = []
    for d in range(len(shape_init)):
        slices.append(slice(center_ext[d]-centre[d],center_ext[d]+shape_init[d]-centre[d]))
    slices = tuple(slices)
    ret = np.zeros(shape_ext)
    ret[slices] += val
    return ret

def interpolate(val,shape_int):
    """
    Interpolate field to desired shape.
    """
    from scipy.interpolate import RegularGridInterpolator
    xyz = tuple(np.arange(0,sh) for sh in val.shape)
    xyx_points = (np.linspace(_x[0],_x[-1],sh) for _x,sh in zip(xyz,shape_int))
    points = np.vstack([X.ravel() for X in np.meshgrid(*xyx_points,indexing='ij')]).T
    ret = RegularGridInterpolator(xyz,val)(points).reshape(shape_int)
    return ret

def apodize(val,radius,slope,dL):
    """ Apodize a cluster about the grid centre. """
    xyz = tuple(np.arange(-sh//2,sh//2)*dL for sh in val.shape)
    XYZ = np.meshgrid(*xyz,indexing='ij')
    R = np.sum([X**2 for X in XYZ],axis=0)**.5
    ret = val * 1/(1+np.exp(slope*(R-radius)))
    return ret

def crop_domain(val,shape_cr):
    """ Crop a domain about the grid centre. """
    slices = []
    for sh,sh_cr in zip(val.shape,shape_cr):
        cr = (sh-sh_cr)//2
        slices.append(slice(cr,sh_cr+cr))
    slices = tuple(slices)
    val = val[slices]
    return val[slices]

def generate_vals(varNames):
    """ Generator for ENZO itasca IA2 fields for the inital snapshot. """
    s = snapshots[0]
    for varName in varNames:
        if varName=='rho':
            val = ia2.get_val('Density',s,c=c) / rho_0
        elif varName=='prs':
            rho = ia2.get_val('Density',s,c=c) / rho_0
            temp = ia2.get_val('Temperature',s,c=c) / K_0
            val = rho*temp/mu #TODO is this the right conversion?
        elif varName=='dm':
            val = ia2.get_val('Dark_Matter_Density',s,c=c) / rho_0
        elif varName=='vx1':
            val = ia2.get_val('x-velocity',s,c=c) / v_0
        elif varName=='vx2':
            val = ia2.get_val('y-velocity',s,c=c) / v_0
        elif varName=='vx3':
            val = ia2.get_val('z-velocity',s,c=c) / v_0
        elif varName=='Bx1':
            val = ia2.get_val('Bx',s,c=c) / B_0
        elif varName=='Bx2':
            val = ia2.get_val('By',s,c=c) / B_0
        elif varName=='Bx3':
            val = ia2.get_val('Bz',s,c=c) / B_0
        else:
            raise ValueError(f"Variable name {varName} is not accounted for.")
        yield val

def generate_dens():
    """ Generator for total density for snapshots. """
    for s in snapshots:
        val = ia2.get_val('Density',s,c=c) + ia2.get_val('Dark_Matter_Density',s,c=c)
        yield val / rho_0

def plot(los):
    """ Plot the final apodised fields. """
    from matplotlib.colors import LogNorm
    from utils.visualise import Style, get_slice
    Style.set()

    plotNames = ['vx1','vx2','vx3','Bx1','Bx2','Bx3','rho','prs']
    for i in range(len(snapshots)):
        plotNames.append(f'pot{i}')
    def load_vals():
        for varName in plotNames:
            with open(os.path.join(savePath, f'{varName}0.dbl'),'rb') as f_o:
                val = np.array(struct.unpack('<'+'d'*n_points,f_o.read())).reshape(shape_final).T
            yield val
    vals = load_vals()

    nrows = 3*np.ceil(len(plotNames)/3)
    fig = plt.figure(figsize=(10,nrows),tight_layout=True)
    for i,(varName,val) in enumerate(zip(plotNames,vals)):
        varkwargs = Style.get_varkwargs(varName)
        ax = fig.add_subplot(nrows,3,i+1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(varkwargs["title"])
        val = get_slice(val,slc=val.shape[0]//2,los=los,c=1)
        norm = None
        if varName in {'rho','dm','prs'}: #TODO
            norm = LogNorm()
        if varName.startswith('vx') or varName.startswith('Bx'):
            vmax = np.abs(val).max()
            vmin = -vmax
        else: vmin, vmax = None, None

        cmap = varkwargs["cmap"]
        im = ax.pcolormesh(val,cmap=cmap,vmin=vmin,vmax=vmax,norm=norm)
        cbar = Style.get_cbar(im,fig,ax,size=.1)
        # cbar.ax.set_title()

    plt.savefig(os.path.join(savePath,f'cluster_{los}.png'),dpi=300)
    plt.show()

if __name__=="__main__":
    
    import os, struct

    from utils.ia2 import IA2, IA2Data
    from utils.clusters import calc_avgvelocity
    from utils.gravity import calc_gravpot

    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("cluster", help="Which IA2 cluster to process.",type=str)
    parser.add_argument("snapshots", help="Which snapshots to animate.")
    parser.add_argument("savePath", help="Where to save ICs.",type=str)
    args = parser.parse_args()
    cluster = args.cluster
    snapshots = [int(s) for s in args.snapshots.removeprefix('[').removesuffix(']').split(',')]
    savePath = args.savePath

    # Load interpolation arguments from topluto.info
    read_topluto(savePath)
    # Load cluster.
    ia2 = IA2Data(cluster)
    # Load cluster units.
    define_units(ia2,Delta)

    # Load time varying cluster properties
    centres = tuple(ia2.get_centre(s,c=c) for s in snapshots)
    redshifts = tuple(IA2.get_redshift(s) for s in snapshots)
    dLs = tuple(IA2.dL * c / (1+z) / L_0 for z in redshifts)

    # Find the smallest volume that will excapsulate the cluster centre at all snapshots.
    L_max = []
    for s,z,cen in zip(snapshots,redshifts,centres):
        L = ia2.dim*IA2.dL / (1+z) / L_0
        cen = tuple(_cen*c*IA2.dL / (1+z) / L_0 for _cen in cen)
        L_max.append(max(tuple(max(_cen,L-_cen) for _cen in cen)) * 2)
    L_max = max(L_max) # Cluster units.

    # Generate raw field values.
    varNames = ["prs","Bx1","Bx2","Bx3"]
    rhovxNames = ["rho","vx1","vx2","vx3"]
    vals_init = generate_vals(varNames)
    rhovxs_init = generate_vals(varNames)
    dens_init = generate_dens()
    # Centre and extend domains.
    shapes_ext = tuple((int(L_max*(1+z)/(IA2.dL/L_0))//c,)*3 for z in redshifts)
    vals_ext = (extend_domain(val,shapes_ext[0],centres[0]) for val in vals_init)
    rhovxs_ext = (extend_domain(val,shapes_ext[0],centres[0]) for val in rhovxs_init)
    dens_ext = (extend_domain(val,shapes_ext[i],centres[i]) for i,val in enumerate(dens_init))
    # Interpolate them to the same dimensions.
    shape_inp = shapes_ext[-1]
    vals_inp = (interpolate(val,shape_inp) for val in vals_ext)
    rhovxs_inp = (interpolate(val,shape_inp) for val in rhovxs_ext)
    dens_inp = (interpolate(val,shape_inp) for val in dens_ext)
    dL = L_max/shape_inp[0] # Cluster units.
    
    # Remove average velocity.
    rho_inp = next(rhovxs_inp)
    vxs_sh = (vx-calc_avgvelocity(rho_inp,vx) for vx in rhovxs_inp)
    rhovxs_sh = (rho_inp if Name=="rho" else next(vxs_sh) for Name in rhovxNames)
    # Calculate gravitational potentials.
    pot = (calc_gravpot(dens,(dL,)*3,G=G_0) for dens in dens_inp)

    # Apodize all fields
    vals_apd = (apodize(val,apd_radius,apd_slope,dL) for val in vals_inp)
    vxs_apd = (apodize(val,apd_radius,apd_slope,dL) for val in vxs_sh)
    # Crop all fields
    shape_cr = (int(L_final//dL),)*3
    vals_cr = (crop_domain(val,shape_cr) for val in vals_apd)
    vxs_cr = (crop_domain(val,shape_cr) for val in vxs_apd)
    pot_cr = (crop_domain(val,shape_cr) for val in pot)
    # Interpolate to final dimension
    vals_final = (interpolate(val,shape_final) for val in vals_cr)
    vxs_final = (interpolate(val,shape_final) for val in vxs_cr)
    pot_final = (interpolate(val,shape_final) for val in pot_cr)

    # Save final fields!
    n_points = np.prod(shape_final)
    for varName,val in zip(varNames,vals_final):
        with open(os.path.join(savePath, f'{varName}0.dbl'),'wb') as f_o:
            f_o.write(struct.pack('<'+'d'*n_points,*(val.T.flatten())))
    for varName,val in zip(rhovxNames[1:],vxs_final):
        with open(os.path.join(savePath, f'{varName}0.dbl'),'wb') as f_o:
            f_o.write(struct.pack('<'+'d'*n_points,*(val.T.flatten())))
    # Save potentials velocities
    for i,val in enumerate(pot_final):
        val -= val.max()
        with open(os.path.join(savePath, f'pot{i}0.dbl'),'wb') as f_o:
            f_o.write(struct.pack('<'+'d'*n_points,*(val.T.flatten())))
    # Save grid.
    DIM = 3
    L = (shape_cr[0]*dL,)*3
    with open(os.path.join(savePath, 'grid0.out'),'w') as f_o:
        f_o.write("# GEOMETRY:   CARTESIAN\n")
        for d in range(DIM):
            f_o.write(f"{shape_final[d]}\n")
            for i in range(shape_final[d]):
                xL = -0.5 + (i-0.5)/(shape_final[d] - 1.)
                xR = -0.5 + (i+0.5)/(shape_final[d] - 1.)
                xL *= L[d]
                xR *= L[d]
                f_o.write("{:d}   {:12.6e}  {:12.6e}\n".format(i+1, xL, xR))

    # Save information.
    H2_a = lambda a : IA2.H_0**2 * (IA2.Omega_r0/a**4 + (IA2.Omega_bm0+IA2.Omega_dm0)/a**3 + IA2.Omega_Lambda0 + IA2.Omega_k0/a**2)
    H2_z = lambda z : H2_a(1/(1+z))
    rho_cr_z = lambda z: 3*H2_z(z)/(8*np.pi*CGS.G)
    rho_min = rho_cr_z(IA2.get_redshift(snapshots[0]))
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
    with open(os.path.join(savePath,'info.txt'),'w') as info:
        info.write(f"L = {np.array(L)/L_0}, L/2 = {np.array(L)/L_0/2}\n")
        info.write(f"rho_0 = {rho_0}\n")
        info.write(f"L_0 = {L_0}\n")
        info.write(f"v_0 = {v_0}\n")
        info.write(f"shape = {shape_final}\n")
        info.write("\n")
        # info.write(f"rho_cr = {rho_min/rho_0}\n")
        # info.write(f"prs_cr = {prs_min/p_0}\n")
        info.write("\n")
        info.write(f"times = {times/t_0}\n")

    plot('x')
    plot('y')
    plot('z')