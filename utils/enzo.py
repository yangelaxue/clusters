"""
Useful data to deal with ENZO and IA2 data.
Rewritten from a February 2026 file.

Author: Angela Xue
Date: March 2026
"""

import numpy as np
import h5py, os

class IA2:

    cd = 2.66e-30
    cv = 6.47e9
    cb = (cd*4*np.pi)**.5 * cv
    dL = 3.95 # comoving kpc
    dim = 1280
    gamma = 5/3

    fields = (
        'Density',
        'Dark_Matter_Density',
        'Temperature',
        'x-velocity',
        'y-velocity',
        'z-velocity',
        'Bx',
        'By',
        'Bz',
    )

    snapshots = (3,4,5,6,7,8,9,10,11,12,13,14,15)
    redshifts = (2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.02,0.01)
    hdf5Name = '8mdmd001R_{}_{}'

class IA2Data:

    def __init__(self,Path):
        self.Path = Path

    def get_val(self,varName,snapshot=None,redshift=None,slc=-1,los=0,c=1,units='cgs'):
        """
        Load field data of a single field from a single snapshot.
        Return it in physical CGS units.
        """

        _const = 1 # Constant to return field values in physical CGS units.
        hdf5Name = os.path.join(self.Path,IA2.hdf5Name)

        # Determine the snapshot to load.
        assert snapshot or redshift, "At least one of snapshot or redshift must be provided."
        if not snapshot:
            snapshot = get_snapshot(redshift)
        z = get_redshift(snapshot)
        snapshot_str = ('00'+str(snapshot))[-3:]
        hdf5Name = hdf5Name.format('{}',snapshot_str)
        
        if varName in {'Density','Dark_Matter_Density','Temperature'}:
            hdf5Name = hdf5Name.format('dt')
        elif varName in {'x-velocity','y-velocity','z-velocity'}:
            hdf5Name = hdf5Name.format('v')
        elif varName in {'Bx','By','Bz'}:
            hdf5Name = hdf5Name.format('b')
        else:
            raise ValueError(f'Variable {varName} is not recognised.')
        
        if units.lower()=='cgs':
            if 'Density' in varName:
                _const *= IA2.cd*(1+z)**3
            elif varName in {'x-velocity','y-velocity','z-velocity'}:
                _const *= IA2.cv
            elif varName in {'Bx','By','Bz'}:
                _const *= IA2.cb*(1+z)**2
        elif units.lower()=='code':
            _const = 1
        else:
            raise ValueError(f"Units {units} is not supported.")

        ret = read_hdf5(hdf5Name,varName,slc,los,c) * _const
        return ret
    
    def gen_vals(self,varName,snapshots=None,redshifts=None,**kwargs):
        """ See get_val """

        if snapshots==None:
            assert not redshifts==None, "Must provide list of snapshots or redshifts"
            snapshots = [get_snapshot(z) for z in redshifts]
        elif snapshots=='all':
            snapshots = list(IA2.snapshots)
        elif type(snapshots)==int:
            snapshots = [snapshots]
        else:
            assert type(snapshots)==list, "Given snapshots must be an integer, list or 'all'."

        if type(slcs)!=list:
            slcs = [slcs,] * len(snapshots)
        else:
            assert type(slcs)==list, "Must ensure that a list of slices is passed."
        
        return (self.get_val(varName,snapshot=s,**kwargs) for s,slc in zip(snapshots,slcs))

    def get_centre(self,snapshot=None,redshift=None,c=1):

        dm = self.get_val('Dark_Matter_Density',snapshot,redshift,c=c)
        rho = self.get_val('Density',snapshot,redshift,c=c)

        dens = dm + rho

        centre_slice = np.where(dens==np.nanmax(dens))
        centre = centre_slice[0][0], centre_slice[1][0], centre_slice[2][0]

        return centre
    
    def save_derivedfield(self,snapshot=None,redshift=None,c=1):
        """
        Save derived field to minimise repeated calculations.

        Implemented: (twice) specific energy.
        """

        vx = (self.get_val(field,snapshot,redshift,c=c) for field in ('x-velocity','y-velocity','z-velocity'))
        v2 = np.zeros((IA2.dim//c,)*3)
        for _vx in vx:
            v2 += _vx**2

        with h5py.File(os.path.join(self.Path,'derivedfields.h5'), 'a') as f:
            f.create_dataset(f"{snapshot}/twice-specific-energy", data=v2)
    
def read_hdf5(hdf5Name,varName,slc=-1,los=0,c=1):
    """
    Loads H5 data from a given file and speficied variable.
    Such variables include:
        'Density' 'Dark_Matter_Density' 'Temperature' 'x-velocity' 'y-velocity' 'z-velocity' 'Bx' 'By' 'Bz'
    
    Parameters
    ----------
        hdf5Name : str
            File name (and location) of the .h5 file to load.
        varName : str
            Name of the variable which to be returned.
        slc : int
            If slc>-1, return the specified slice.
        los : int, str
            If only a slice is to be returned, determine the axis of the slice.
        c : int
            Return a subset of data values at even intervals of c, to reduce space needed.

    
    Returns
    -------
        ret : np.ndarray
            List of the grid of variables returned in the same order as var.
    """

    assert varName in IA2.fields, f"Variable varName must be in {IA2.fields}"

    with h5py.File(hdf5Name, 'r') as f:
        if slc>-1:
            if los in (0,'x','X'):
                ret = f[varName][::c,::c,::c][slc]
            elif los in (0,'y','Y'):
                ret = f[varName][::c,::c,::c][:,slc]
            elif los in (0,'z','Z'):
                ret = f[varName][::c,::c,::c][:,:,slc]
            else:
                ValueError("Axis must be 0, x, X, 1, y, Y, 2, z, or Z.")
        else:
            ret = f[varName][::c,::c,::c]

    return ret

def get_snapshot(z):
    ddict = {_z : _s for _z,_s in zip(IA2.redshifts,IA2.snapshots)}
    return ddict[z]

def get_redshift(s):
    ddict = {_s : _z for _z,_s in zip(IA2.redshifts,IA2.snapshots)}
    return ddict[s]

if __name__=="__main__":

    Path = "/media/yange/MyDrive/2024PhDData/EnzoIA2/E14-R"

    ia2 = IA2Data(Path)
    ia2.save_derivedfield(snapshot=15,c=5)
    # for s in IA2.snapshots:
    #     try:
    #         ia2.save_derivedfield(snapshot=s,c=5)
    #     except:
    #         print(f"Snapshot {s} not found. Continue.")