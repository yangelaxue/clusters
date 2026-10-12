"""
Classes and functions that deal with the ENZO Itasca IA2 clusters.

Author: Angela Xue
Date: Feb 2026
Revised: Oct 2026
"""

import numpy as np
import os, h5py
from utils.units import CGS

class IA2:
    """ Information for all clusters. """

    # Cosmological Parameters
    Omega_r0 = 0.
    Omega_bm0 = 0.0441
    Omega_dm0 = 0.2139
    Omega_Lambda0 = 0.742
    Omega_k0 = 0.
    h_0 = 0.72
    H_0 = 100*h_0*1e5/(CGS.pc*1e6) # CGS units.
    dL = 3.95 * (1e3*CGS.pc) # Comoving cell width.

    # Conversion to CGS units
    cd = 2.66e-30
    cv = 6.47e9
    cb = (cd*4*np.pi)**.5 * cv

    # Fields stored.
    fields = {
        "Density",
        "Dark_Matter_Density",
        "Temperature",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Bx",
        "By",
        "Bz",
    }
    redshifts = (2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.02,0.01)
    get_snapshot = lambda z : {_z : _s for _z,_s in zip(IA2.redshifts,IA2.snapshots["all"])}[z]
    get_redshift = lambda s : {_s : _z for _z,_s in zip(IA2.redshifts,IA2.snapshots["all"])}[s]
    # Cluster specific snapshots
    snapshots = {
        "all" : (3,4,5,6,7,8,9,10,11,12,13,14,15),
        "E14" : (3,4,5,6,7,8,10,11,12,13,14,15),
        "E18B" : (3,4,5,6,7,8,9,10,11,12,13,14,15),
        "E3A" : (3,4,5,6,7,8,9,10,11,12,13,14,15),
        "E5A" : (3,4,5,6,7,8,9,10,11,12,13,14,15),
    }

    #TODO save times.
    hdf5Name = "8mdmd001R_{}_{}"

    # Define path for the clusters
    if os.uname().nodename=="yange":
        Path = "/media/yange/MyDrive/2024PhDData/EnzoIA2/"
    elif "nesi" in os.uname().nodename:
        Path = "/home/xuean/00_nesi_projects/vuw04655_nobackup/IA2"
    else:
        raise SystemError("There is no path for the ENZO Itasca IA2 clusters in this system.")

class IA2Data:
    """ IO of clusters. """

    def __init__(self,cluster):

        self.cluster = cluster
        self.Path = os.path.join(IA2.Path,cluster)

        if cluster=="E3A":
            self.dim = 1024
        else:
            self.dim = 1280
        self.snapshots = IA2.snapshots[cluster]

    def get_val(self,varName,snapshot,units="cgs",**kwargs):
        """
        Load and return field values of a single field from a specified snapshot.
        See get_val for details.

        Will return raw fields or select derived fields.

        Default to CGS units.
        """

        if varName in IA2.fields:
            z = IA2.get_redshift(snapshot)
            _const = 1.
            if units.lower()=="cgs":
                if "Density" in varName:
                    _const *= IA2.cd*(1+z)**3
                elif varName in {"x-velocity","y-velocity","z-velocity"}:
                    _const *= IA2.cv
                elif varName in {"Bx","By","Bz"}:
                    _const *= IA2.cb*(1+z)**2
            elif units.lower().startswith("delta"):
                Delta = int(units.lower().removeprefix("delta"))
                if "Density" in varName:
                    _const *= IA2.cd*(1+z)**3
                    _const /= self.get_units(Delta)[4]
                elif varName in {"x-velocity","y-velocity","z-velocity"}:
                    _const *= IA2.cv
                    _const /= self.get_units(Delta)[2]
                elif varName in {"Bx","By","Bz"}:
                    _const *= IA2.cb*(1+z)**2
            elif units.lower()=="code":
                _const = 1
            else:
                raise ValueError(f"Units {units} is not supported.")
            ret = get_val(self.Path,varName,snapshot,**kwargs) * _const
            return ret
        else:
            assert varName in {"ek","vorticity","B2","curlB"}, f"Derived field {varName} is not recognised."
            if varName=="ek":
                dens = self.get_val("Density",snapshot,units=units,**kwargs)
                dens += self.get_val("Dark_Matter_Density",snapshot,units=units,**kwargs)
                vx = self.get_val("x-velocity",snapshot,units=units,**kwargs)
                vy = self.get_val("y-velocity",snapshot,units=units,**kwargs)
                vz = self.get_val("z-velocity",snapshot,units=units,**kwargs)
                return 0.5 * dens * (vx**2 + vy**2 + vz**2)
            elif varName=="vorticity":
                from utils.gradient import curl_discrete
                vx = self.get_val("x-velocity",snapshot,units=units,**kwargs)
                vy = self.get_val("y-velocity",snapshot,units=units,**kwargs)
                vz = self.get_val("z-velocity",snapshot,units=units,**kwargs)
                return curl_discrete(vx,vy,vz,dxdydz=(IA2.dL,)*3)
            elif varName=="B2":
                Bx = self.get_val("Bx",snapshot,units=units,**kwargs)
                By = self.get_val("By",snapshot,units=units,**kwargs)
                Bz = self.get_val("Bz",snapshot,units=units,**kwargs)
                return Bx**2 + By**2 + Bz**2
            elif varName=="curlB":
                from utils.gradient import curl_discrete
                Bx = self.get_val("Bx",snapshot,units=units,**kwargs)
                By = self.get_val("By",snapshot,units=units,**kwargs)
                Bz = self.get_val("Bz",snapshot,units=units,**kwargs)
                return curl_discrete(Bx,By,Bz,dxdydz=(IA2.dL,)*3)

    def gen_vals(self,varName=None,snapshot=None,varNames=None,snapshots=None,units="cgs",**kwargs):
        """
        Generator for field values.
            For a single field for many snapshots,
        or
            For a many fields for a single snapshot.
        """
        assert varNames or snapshots, "Only generate fields for a single snapshot, or snapshots for a single field."
        if varNames:
            for varName in varNames:
                yield self.get_val(varName,snapshot,units,**kwargs)
        if snapshots:
            for s in snapshots:
                yield self.get_val(varName,s,units,**kwargs)

    def get_centre(self,snapshot,c=1):
        """ Return the center of mass index saved. """
        with open(os.path.join(self.Path,"centreindices.txt"),"r") as f:
            for line in f.readlines():
                line = [_line for _line in line.split(" ") if _line]
                s, idx, idy, idz = line
                if s=="s":
                    continue
                elif int(s)==snapshot:
                    centreidx = (int(float(idx)//c), int(float(idy)//c), int(float(idz)//c))
                    return centreidx
        print(f"Could not locate centre index for snapshot={snapshot}.")
        return None

    def get_units(self,Delta):
        """
        Return cluster units in CGS units.
        Units returned are: r_Delta, M_Delta, c_s, t_s, centre density.
        """
        with open(os.path.join(self.Path,"clusterunits.txt"),"r") as f:
            for line in f.readlines():
                line = [_line for _line in line.split(" ") if _line]
                _Delta, r_Delta, M_Delta, c_s, t_s = line
                if _Delta=="Delta":
                    continue
                if int(_Delta)==Delta:
                    r_Delta, M_Delta, c_s, t_s = float(r_Delta), float(M_Delta), float(c_s), float(t_s)
                    break
        with open(os.path.join(self.Path,"centredensity.txt"),"r") as f:
            for line in f.readlines():
                line = [_line for _line in line.split(" ") if _line]
                s, density = line
                if s=="s":
                    continue
                if int(s)==15:
                    density = float(density)
                    return r_Delta, M_Delta, c_s, t_s, density
        print(f"Could not locate cluster parameters or centre density.")
        return None

    def get_fnames(self):
        """
        Returns list of data files in order of snapshots.
        """

        varstr = ["dt","v","b"]

        fNames = {_varstr : [] for _varstr in varstr}
        for fName in os.listdir(self.Path):
            for _varstr in varstr:
                if fName.startswith(f"8mdmd001R_{_varstr}"):
                    fNames[_varstr].append(fName)
        for _varstr in varstr:
            idxes = np.argsort([int(fName[-3:]) for fName in fNames[_varstr]])
            fNames.update({_varstr : [fNames[_varstr][idx] for idx in idxes]})

        return fNames

def get_val(Path,varName,snapshot,c=1):
    """
    Loads H5 data from a given file and speficied variable.
    Such variables include:
        'Density' 'Dark_Matter_Density' 'Temperature' 'x-velocity' 'y-velocity' 'z-velocity' 'Bx' 'By' 'Bz'
    
    Parameters
    ----------
        Path : str
            Location of the stored field values.
        varName : str
            Name of the variable which to be returned.
        snapshot : int
            Snapshot which to load the field values for.
        c : int
            Return a subset of data values at even intervals of c, to reduce space needed.
    Returns
    -------
        ret : np.ndarray
            Field values (code units).
    """

    snapshot_str = ("00"+str(snapshot))[-3:]
    hdf5Name = os.path.join(Path,IA2.hdf5Name).format("{}",snapshot_str)
    if varName in {"Density","Dark_Matter_Density","Temperature"}:
        hdf5Name = hdf5Name.format("dt")
    elif varName in {"x-velocity","y-velocity","z-velocity"}:
        hdf5Name = hdf5Name.format("v")
    elif varName in {"Bx","By","Bz"}:
        hdf5Name = hdf5Name.format("b")
    else:
        raise ValueError(f"Variable {varName} is not recognised.")

    with h5py.File(hdf5Name,"r") as f:
        ret = f[varName][::c,::c,::c]
    return ret