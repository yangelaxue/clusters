"""
This script will process IA2 clusters to be PLUTO initial conditions automatically from start to finish.

Functions include:
    define_units - defines code units to use in PLUTO.
    load_rawinit, load_rawrhos, load_rawdms : generators for field values.

Author: Angela Xue
Data: June 2026
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def define_units(Delta):
    """ Define code units to use in PLUTO. These are cluster units. """
    params_fName = os.path.join(ia2.Path,f'cluster_{Delta}.txt')
    r_Delta, M_Delta, rho_Delta, c_s, t_s = np.loadtxt(params_fName,skiprows=1,max_rows=1,delimiter=',')
    rho_cs = np.loadtxt(os.path.join(ia2.Path,f'rho_c.txt'),skiprows=1)
    idx, = np.where(rho_cs==15)[0]
    rho_c = rho_cs[idx,1]

    global L_0, v_0, rho_0
    global t_0, B_0, p_0
    global K_0, G_0

    L_0 = r_Delta
    v_0 = c_s
    rho_0 = rho_c
    
    t_0 = L_0/v_0 # time unit in seconds
    B_0 = v_0 * (4*np.pi*rho_0)**.5 # Magnetic field unit in Gauss (G)
    p_0 = rho_0*v_0**2 # Pressure unit in dyne/cm^2
    
    K_0 = CGS.mp * v_0**2 / CGS.kB
    G_0 = CGS.G / (1/rho_0/t_0**2)

def load_rawinit():
    """ Generator for initial field values, sans velcity fields. """

    s = snapshots[0]
    for varName in varNames:
        if varName=='rho':
            val = ia2.get_val('Density',s,c=c) / rho_0
        elif varName=='prs':
            rho = ia2.get_val('Density',s,c=c) / rho_0
            temp = ia2.get_val('Temperature',s,c=c) / K_0
            val = rho*temp/mu
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

def load_rawvxs():
    """ Generator for initial density and velocity field values. """
    
    s = snapshots[0]
    for varName in vxNames:
        if varName=='rho':
            val = ia2.get_val('Density',s,c=c) / rho_0
        elif varName=='vx1':
            val = ia2.get_val('x-velocity',s,c=c) / v_0
        elif varName=='vx2':
            val = ia2.get_val('y-velocity',s,c=c) / v_0
        elif varName=='vx3':
            val = ia2.get_val('z-velocity',s,c=c) / v_0
        else:
            raise ValueError(f"Variable name {varName} is not accounted for.")
        yield val

def load_rawrhos():
    """ Generator for density values. """
    for s in snapshots:
        yield ia2.get_val("Density",snapshot=s,c=c,units='cgs') / rho_0

def load_rawdms():
    """ Generator for dark matter density values. """
    for s in snapshots:
        yield ia2.get_val("Dark_Matter_Density",snapshot=s,c=c,units='cgs') / rho_0

def centre_domain(centre, val):
    """
    Centre field values on center of mass by extending domain.
    """

    shape = val.shape
    pads = []
    for cen,sh in zip(centre,shape):
        pads.append([])
        if cen<sh//2:
            pads[-1].append(sh-cen*2)
            pads[-1].append(0)
        else:
            pads[-1].append(0)
            pads[-1].append(cen*2-sh)

    shape_final = max([sh+max(pad) for pad,sh in zip(pads,val.shape)])
    
    ret = np.zeros((shape_final,)*val.ndim)
    slices = []
    for pad,sh in zip(pads,val.shape):
        _pad = (shape_final - (sh + sum(pad)))//2
        pad[0] += _pad
        pad[1] += _pad
        slices.append(slice(pad[0],shape_final-pad[1]))
    slices = tuple(slices)
    ret[slices] = val
    
    return ret

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

def interpolate(shape_int,val):
    """
    Interpolate field to desired shape.
    """

    xyz = tuple(np.arange(0,sh) for sh in val.shape)

    xyx_points = (np.linspace(_x[0],_x[-1],sh) for _x,sh in zip(xyz,shape_int))
    points = np.vstack([X.ravel() for X in np.meshgrid(*xyx_points,indexing='ij')]).T
    ret = RegularGridInterpolator(xyz,val)(points).reshape(shape_int)

    return ret

def get_R(shape,dL):
    xyz = tuple(np.arange(-sh//2,sh//2)*dL for sh in shape)
    XYZ = np.meshgrid(*xyz,indexing='ij')
    R = np.sum([X**2 for X in XYZ],axis=0)**.5
    return R

def sigmoid(x,radius,slope):
    return 1/(1+np.exp(slope*(x-radius)))

def crop_domain(shape_cr,val):

    slices = []
    for sh_cr,sh in zip(shape_cr,val.shape):
        cr = (sh-sh_cr)//2
        slices.append(slice(cr,sh_cr+cr))
    slices = tuple(slices)

    return val[slices]

def get_times():
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmo = FlatLambdaCDM(H0=IA2.h_0*100*u.km/u.s/u.Mpc, Om0=IA2.Omega_dm0+IA2.Omega_bm0)

    times = []
    for s in snapshots:
        z = get_redshift(s)
        times.append(cosmo.age(z).value*(1e9*CGS.yr))
    times = np.array(times)
    return times-times.min()

def plot():

    plotNames = ['vx1','vx2','vx3','Bx1','Bx2','Bx3','rho','prs']
    for i in range(len(snapshots)):
        plotNames.append(f'pot{i}')
    def load_vals():
        for varName in plotNames:
            with open(os.path.join(savePath, f'{varName}0.dbl'),'rb') as f_o:
                val = np.array(struct.unpack('<'+'d'*n_points,f_o.read())).reshape(shape_cr).T
            yield val
    vals = load_vals()

    fig = plt.figure(figsize=(10,3*np.ceil(len(plotNames)/3)),tight_layout=True)

    for i,(varName,val) in enumerate(zip(plotNames,vals)):

        ax = fig.add_subplot(4,3,i+1)
        ax.set_aspect(1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(varName)
        val = val[shape_cr[0]//2]
        norm = None
        if varName in {'rho','dm','prs'}: #TODO
            norm = LogNorm()
        if varName.startswith('vx') or varName.startswith('Bx'):
            vmax = np.abs(val).max()
            vmin = -vmax
        else: vmin, vmax = None, None

        cmap = DefaultStyle.get_varkwargs(varName)['cbar_cmap']
        im = ax.pcolormesh(val,cmap=cmap,vmin=vmin,vmax=vmax,norm=norm)
        cbar = get_cbar(im,fig,ax,cbar_pad=0,cbar_size=.1)

    plt.savefig(os.path.join(savePath,'cluster.png'),dpi=300)

def plot_apodize():
    x, y = (np.arange(-shape_int[0]//2,shape_int[0]//2)*L_max/L_0/shape_int[0],)*2
    fig, ax = plt.subplots(1,1,figsize=(5,4),tight_layout=True)
    im = ax.pcolormesh(x,y,apodize[shape_int[0]//2],cmap='Greys',vmin=0.,vmax=1.)
    ax.set_aspect(1)
    cbar = get_cbar(im,fig,ax)
    plt.savefig(os.path.join(savePath,'apodize.png'),dpi=300)

H2_a = lambda a : IA2.H_0**2 * (IA2.Omega_r0/a**4 + (IA2.Omega_bm0+IA2.Omega_dm0)/a**3 + IA2.Omega_Lambda0 + IA2.Omega_k0/a**2)
H2_z = lambda z : H2_a(1/(1+z))
rho_cr_z = lambda z: 3*H2_z(z)/(8*np.pi*CGS.G)

if __name__=="__main__":
    
    import argparse, os, struct
    from utils.enzo import IA2, IA2Data, get_redshift
    from utils.units import CGS
    from utils.gravity import get_gravpot
    from utils.clusters import calc_avgvelocity
    from utils.visualize import DefaultStyle, get_cbar, set_style
    set_style()

    parser = argparse.ArgumentParser("main")
    parser.add_argument("cluster", help="Which IA2 cluster to process.", type=str)
    parser.add_argument("savePath", help="Where to save ICs.", type=str)
    parser.add_argument("snapshots", help="Which snapshots to animate.")

    #%% Parse arguments

    args = parser.parse_args()
    cluster = args.cluster
    savePath = args.savePath
    if not os.path.exists(savePath):
        os.makedirs(savePath)
    snapshots = args.snapshots if args.snapshots else 'all'
    if snapshots=='all':
        snapshots = list(IA2.snapshots)
    else:
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]

    #%% Load cluster class and information.

    IA2Paths = {
        "E14" : os.path.join(os.environ['IA2_DIR'],"E14"),
        "E18B" : os.path.join(os.environ['IA2_DIR'],"E18B"),
        "E3A" : os.path.join(os.environ['IA2_DIR'],"E3A"),
        "E5A" : os.path.join(os.environ['IA2_DIR'],"E5A"),
    }
    ia2 = IA2Data(IA2Paths[cluster])

    # Load parameters from file. Else, default.
    c, Delta, mu, radius, slope, shape_cr = 40, 200, 0.5, 0.5, 16., ()
    with open(os.path.join(savePath,'topluto.info'),'r') as info:
        for line in info.readlines():
            line = line.split('#')[0].replace(" ","").replace("\n","").replace("\t","")
            if line.startswith('c'):
                c = int(line.split('=')[-1])
                continue
            if line.startswith('Delta'):
                Delta = int(line.split('=')[-1])
                continue
            if line.startswith('mu'):
                mu = float(line.split('=')[-1])
                continue
            if line.startswith('radius'):
                radius = float(line.split('=')[-1])
                continue
            if line.startswith('slope'):
                slope = float(line.split('=')[-1])
                continue
            if line.startswith('shape_cr'):
                line = line.split('=')[-1]
                try:
                    shape_cr = (int(line),)*3
                except:
                    shape_cr = tuple(int(sh) for sh in line.removeprefix('(').removesuffix(')').split(','))
                continue
    shape_cr = shape_cr if shape_cr else (ia2.dim//c,)*3

    # Load important.
    centres = tuple(ia2.get_centre(s) for s in snapshots)
    redshifts = tuple(get_redshift(s) for s in snapshots)
    dLs = tuple(IA2.dL * (1e3*CGS.pc) * c / (1+z) for z in redshifts)
    define_units(Delta)

    #%% Load field values.
    varNames = ['rho','prs','Bx1','Bx2','Bx3']
    vxNames = ['rho','vx1','vx2','vx3']
    vals_raw = load_rawinit()
    vxs_raw = load_rawvxs()
    rawrhos = load_rawrhos()
    rawdms = load_rawdms()

    #%% Centre fields.
    vals_cen = (centre_domain(tuple(cen//c for cen in centres[0]),val) for val in vals_raw)
    vxs_cen = (centre_domain(tuple(cen//c for cen in centres[0]),val) for val in vxs_raw)
    rhos_cen = (centre_domain(tuple(cen//c for cen in centre),val) for val,centre in zip(rawrhos,centres))
    dms_cen = (centre_domain(tuple(cen//c for cen in centre),val) for val,centre in zip(rawdms,centres))

    #%% Calculate maximum domain length before final crop.
    Ls = []
    for i,s in enumerate(snapshots):
        Ls.append(max([max(cen//c,ia2.dim//c-cen//c)*2 for cen in centres[i]]) * dLs[i])
    L_max = max(Ls)

    #%% Extend all fields to same physical length.
    vals_ext = (extend_domain((int(L_max/Ls[0]*val.shape[0]),)*3,val) for i,val in enumerate(vals_cen))
    vxs_ext = (extend_domain((int(L_max/Ls[0]*val.shape[0]),)*3,val) for i,val in enumerate(vxs_cen))
    rhos_ext = (extend_domain((int(L_max/Ls[i]*val.shape[0]),)*3,val) for i,val in enumerate(rhos_cen))
    dms_ext = (extend_domain((int(L_max/Ls[i]*val.shape[0]),)*3,val) for i,val in enumerate(dms_cen))

    #%% Interpolate all grids to same dimension.
    shape_int = []
    for i,s in enumerate(snapshots):
        shape_int.append(round(max([max(cen//c,ia2.dim//c-cen//c)*2 for cen in centres[i]]) * L_max/Ls[i]))
    shape_int = (max(shape_int),)*3
    
    vals_int = (interpolate(shape_int,val) for val in vals_ext)
    vxs_int = (interpolate(shape_int,val) for val in vxs_ext)
    rhos_int = (interpolate(shape_int,val) for val in rhos_ext)
    dms_int = (interpolate(shape_int,val) for val in dms_ext)

    #%% Calculate dL.
    dL = L_max/shape_int[0]

    #%% Apodize fields.
    R = get_R(shape_int,dL/L_0)
    apodize = sigmoid(R,radius,slope)
    vals_apd = (apodize*val for val in vals_int)
    vxs_apd = (apodize*val for val in vxs_int)
    rhos_apd = (apodize*val for val in rhos_int)
    dms_apd = (apodize*val for val in dms_int)

    plot_apodize()

    #%% Calculate gravitational potentials.
    slices_int = tuple(slice(sh//2,sh*3//2) for sh in shape_int)
    gravpots = (get_gravpot(extend_domain(tuple(sh*2 for sh in shape_int),rho+dm),(dL/L_0,)*3,G=G_0)[slices_int] for rho,dm in zip(rhos_apd,dms_apd))

    #%% Crop final fields.
    vals_cr = (crop_domain(shape_cr,val) for val in vals_apd)
    vxs_cr = (crop_domain(shape_cr,val) for val in vxs_apd)
    gravpots_cr = (crop_domain(shape_cr,val)-val.max() for val in gravpots)

    #%% Save fields. Final loop through all values.
    n_points = np.prod(shape_cr)

    # Save fields besides (non) shifted velocity.
    for varName,val in zip(varNames,vals_cr):
        if varName=='prs':
            prs_min = val[val>0].min()*p_0 # Save minimum pressure value.
        with open(os.path.join(savePath, f'{varName}0.dbl'),'wb') as f_o:   
            f_o.write(struct.pack('<'+'d'*n_points,*(val.T.flatten())))
    
    # Save shifted velocities
    rho = next(vxs_cr)
    vxs_sh = (vx-calc_avgvelocity(rho,vx) for vx in vxs_cr)
    for i,vx in enumerate(vxs_sh):
        with open(os.path.join(savePath, f'vx{i+1}0.dbl'),'wb') as f_o:   
            f_o.write(struct.pack('<'+'d'*n_points,*(vx.T.flatten())))
    # Save gravitational potentials
    for i,val in enumerate(gravpots_cr):
        with open(os.path.join(savePath, f'pot{i}0.dbl'),'wb') as f_o:   
            f_o.write(struct.pack('<'+'d'*n_points,*(val.T.flatten())))
    
    #%% Save grid.

    DIM = 3
    L = tuple(dL*sh for sh in shape_cr)
    with open(os.path.join(savePath, 'grid0.out'),'w') as f_o:
        f_o.write("# GEOMETRY:   CARTESIAN\n")
        for d in range(DIM):
            f_o.write(f"{shape_cr[d]}\n")
            for i in range(shape_cr[d]):
                xL = -0.5 + (i-0.5)/(shape_cr[d] - 1.)
                xR = -0.5 + (i+0.5)/(shape_cr[d] - 1.)
                xL *= L[d]
                xR *= L[d]
                f_o.write("{:d}   {:12.6e}  {:12.6e}\n".format(i+1, xL, xR))

    #%% Save information.
    rho_min = rho_cr_z(get_redshift(snapshots[0]))
    times = get_times()
    # prs_min = rho_cr/CGS.mp * CGS.kB * 2.75
    with open(os.path.join(savePath,'info.txt'),'w') as info:
        info.write(f"L = {np.array(L)/L_0}, L/2 = {np.array(L)/L_0/2}\n")
        info.write(f"rho_0 = {rho_0}\n")
        info.write(f"L_0 = {L_0}\n")
        info.write(f"v_0 = {v_0}\n")
        info.write(f"shape = {shape_cr}\n")
        info.write("\n")
        info.write(f"rho_cr = {rho_min/rho_0}\n")
        info.write(f"prs_cr = {prs_min/p_0}\n")
        info.write("\n")
        info.write(f"times = {times/t_0}\n")

    # Load everything and plot.
    plot()