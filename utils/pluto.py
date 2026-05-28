"""
Functions that help with reading and loading data from PLUTO simulations.

Classes
-------
    PlutoData
        Deals with PLUTO data format.
        Loads data based on save Path.

Functions
---------
    get_hdf5Names
        Returns all PLUTO outout file names, which are of the form data.XXXX.dbl.h5
    get_varNames
        Returns all the variables saved by PLUTO run.
    get_times
        Returns all the times corresponding to saved data.
        TODO: I did not double check if it returns slice times or other times...
    get_bbox
        Calculates bounding box of run.
    read_hdf5
        Returns data from specified file with specified variable.
        Can choose which slice or subset of data.
    read_definitions_h
        Return parameter value of speficied parameter.
    read_pluto_ini
        Return parameter value of speficied parameter.

TODO: AMR grid is not accounted for.

Author: Angela Xue
Date: March 2026
"""

import numpy as np
import os, h5py

from utils.units import CGS

class PlutoData:

    def __init__(self,Path:str):

        self.Path = Path
        self._set_metadata()

    def _set_metadata(self,):
        """
        Retrieve the metadata from the file named definitions.h, pluto.ini and derived quantities.
        """

        # Load all data.XXXX.dbl.h5 names and locations.
        self.hdf5Names = [os.path.join(self.Path,hdf5Name) for hdf5Name in get_hdf5Names(self.Path)]

        self.ndim = read_definitions_h(self.Path,'DIMENSIONS')
        self.PHYSICS = read_definitions_h(self.Path,'PHYSICS')
        self.EOS = read_definitions_h(self.Path,'EOS')
        if self.EOS=='IDEAL':
            self.gamma = 5/3
        try:
            self.FORCED_TURB = read_definitions_h(self.Path,'FORCED_TURB')
        except:
            self.FORCED_TURB = 'NO'
        if self.FORCED_TURB=='YES':
            self.FORCED_TURB_DECAY = read_definitions_h(self.Path,'FORCED_TURB_DECAY')
            self.FORCED_TURB_ENERGY = read_definitions_h(self.Path,'FORCED_TURB_ENERGY')
            self.FORCED_TURB_KMIN = read_definitions_h(self.Path,'FORCED_TURB_KMIN')
            self.FORCED_TURB_KMAX = read_definitions_h(self.Path,'FORCED_TURB_KMAX')

        # Units. If not defined, set to default.
        try:
            self.rho_0 = read_definitions_h(self.Path,'UNIT_DENSITY')
        except:
            self.rho_0 = CGS.mp     # proton mass per cm^3
        try:
            self.L_0 = read_definitions_h(self.Path,'UNIT_LENGTH')
        except:
            self.L_0 = CGS.AU       # astronomical unit
        try:
            self.v_0 = read_definitions_h(self.Path,'UNIT_VELOCITY')
        except:
            self.v_0 = 1e5          # 1 km per s
        # Set derived units.
        self.t_0 = self.L_0/self.v_0
        self.p_0 = self.rho_0*self.v_0**2
        self.B_0 = self.v_0*(4*np.pi*self.rho_0)**.5

        # Define grid dimension.
        shape_x, shape_y, shape_z = 1, 1, 1
        shape_x = int(read_pluto_ini(self.Path,'X1-grid')[2])
        shape_y = int(read_pluto_ini(self.Path,'X2-grid')[2])
        shape_z = int(read_pluto_ini(self.Path,'X3-grid')[2])
        self.shape = (shape_x, shape_y, shape_z,)

        self.times = get_times(self.Path)
        self.bbox = get_bbox(self.Path)
        self.varNames = get_varNames(self.Path)

    def print_attr(self,):

        print(f'PHYSICS : {self.PHYSICS}')
        print(f'ndim : {self.ndim}')
        print(f'shape : {self.shape}')
        print(f'bbox : {self.bbox}')
        print(f'_vars : {self._vars}')
        print(f'EOS : {self.EOS}')
        if self.EOS=='IDEAL':
            print(f'gamma : {self.gamma}')
        print()
        print(f'rho_0 : {self.rho_0} {CGS.rho_unit}')
        print(f'L_0 : {self.L_0} {CGS.L_unit}')
        print(f'v_0 : {self.v_0} {CGS.v_unit}')
        print(f't_0 : {self.t_0} {CGS.t_unit}')
        print(f'p_0 : {self.p_0} {CGS.p_unit}')
        print(f'B_0 : {self.B_0} {CGS.B_unit}')
        print()
        print(f'FORCED_TURB : {self.FORCED_TURB}')
        if self.FORCED_TURB=='YES':
            print(f'FORCED_TURB_DECAY : {self.FORCED_TURB_DECAY}')
            print(f'FORCED_TURB_ENERGY : {self.FORCED_TURB_ENERGY}')

    def get_times(self,snapshots='all',units='code'):

        if units=='code':
            _const = 1.
        elif units=='Myr':
            _const = self.t_0/(1e6*CGS.yr)
        elif units=='Gyr':
            _const = self.t_0/(1e9*CGS.yr)
        else:
            print(f"The unit {units} is not supported; defaulting to code units.")
            _const = 1.

        if type(snapshots)==int:
            times = self.times[snapshots]
        elif type(snapshots)==list:
            times = np.array([self.times[i] for i in snapshots])
        elif snapshots=='all':
            times = self.times
        else:
            ValueError(f'snapshots = {snapshots} is not supported.')
        return times*_const
    
    # TODO Move units section outside of class
    def get_x(self,axis:int,units='code'):
        """ Get domain array. """
        
        if axis in {'x','X'}:
            axis = 0
        elif axis in {'y','Y'}:
            axis = 1
        elif axis in {'z','Z'}:
            axis = 2

        x = np.linspace(*self.bbox[axis],self.shape[axis])
        
        _const = 1
        if units=='code':
            return x
        else:
             x *= self.L_0
        if 'pc' in units:
            x /= CGS.pc
        elif 'AU' in units:
            x /= CGS.AU
        else:
            raise NotImplementedError("Only pc and AU lengths are implemented.")
        if units.startswith('k'):
            x /= 1e3
        elif units.startswith('M'):
            x /= 1e6
        elif units.startswith('G'):
            x /= 1e9
        else:
            raise NotImplementedError("Only k, M and G multipliers are implemented.")
        return x
    
    # TODO Move units section outside of class
    def get_dxdydz(self,axis:int,units='code'):
        """ Get domain array. """

        dxdydz = np.array([np.diff(_b)/_s for _b,_s in zip(self.bbox,self.shape)])
        
        _const = 1
        if units=='code':
            return tuple(dx for dx in dxdydz.flatten())
        else:
            dxdydz *= self.L_0
        if 'pc' in units:
            dxdydz /= CGS.pc
        elif 'AU' in units:
            dxdydz /= CGS.AU
        else:
            raise NotImplementedError("Only pc and AU lengths are implemented.")
        if units.startswith('k'):
            dxdydz /= 1e3
        elif units.startswith('M'):
            dxdydz /= 1e6
        elif units.startswith('G'):
            dxdydz /= 1e9
        else:
            raise NotImplementedError("Only k, M and G multipliers are implemented.")
        dxdydz = tuple(dx for dx in dxdydz.flatten())
        return dxdydz

    def get_X(self,X:str):
        """ Get domain mesh. """
        
        if X in {0,'x'}:
            X = 'X'
        elif X in {1,'y'}:
            X = 'Y'
        elif X in {2,'z'}:
            X = 'Z'

        hdf5Name = self.hdf5Names[0]
        return read_hdf5(hdf5Name,'X')
    
    def get_vals(self,varName:str,snapshots='all',units='code',**kwargs):
        """
        Return specified variable at specified snapshot(s).
        """

        snapshots = [snapshots] if type(snapshots)==int else snapshots

        hdf5Names = self.hdf5Names
        if snapshots!='all':
            assert type(snapshots)==list, "snapshots needs to be a list of integers"
            hdf5Names = [self.hdf5Names[s] for s in snapshots]

        if units.lower()=='code':
            const = 1.
        elif units.lower()=='cgs':
            if varName=='rho':
                const = self.rho_0
            elif varName=='prs':
                const = self.p_0
            elif varName.startswith('vx'):
                const = self.v_0
            elif varName.startswith('Bx'):
                const = self.B_0
            else:
                raise NameError(f'Field value with name {varName} does not exist.')
        else:
            raise NotImplementedError("Class can only return field values in code or CGS units.")
        
        return [read_hdf5(hdf5Name,varName,**kwargs)*const for hdf5Name in hdf5Names]
    
    def gen_vals(self,varName,snapshots='all',**kwargs):
        """
        Set a generator for specified variable at specified snapshot(s).

        Parameters
        ----------
            See self.get_vals
        """

        snapshots = [snapshots] if type(snapshots)==int else snapshots
        if snapshots=='all':
            snapshots = [i for i in range(len(self.hdf5Names))]
        
        return (self.get_vals(varName,snapshot,**kwargs)[0] for snapshot in snapshots)

# Functions that can be used outside the class.

def get_hdf5Names(Path:str):
    """
    Load all HDF5 file names in a directory in ascending order.
    A maximum of 1000 output files are produced by PLUTO.
    """

    f_names = [f_name for f_name in os.listdir(Path) if f_name.endswith('.dbl.h5')]
    f_nums = np.sort([int(f_name[5:9]) for f_name in f_names])
    sorted_hdf5Names = [f'data.{('000'+str(f_num))[-4:]}.dbl.h5' for f_num in f_nums]

    return sorted_hdf5Names

def get_varNames(Path:str):
    """ Retrieve field variable names. """

    fName = get_hdf5Names(Path)[0]

    with h5py.File(os.path.join(Path,fName)) as f:
        varNames = list(f['Timestep_0']['vars'].keys())
    return varNames

def get_times(Path:str):
    """ Get all times from the file 'dbl.h5.out'. """

    return np.loadtxt(os.path.join(Path,'dbl.h5.out'), usecols=1)

def get_bbox(Path:str):
    """ Retrive bounding box. """
    bbox = []
    with open(os.path.join(Path,'grid.out')) as f:
        for line in f.readlines()[:20]: # only need to parse the top of the file.
            if line.startswith('# X1: ') or line.startswith('# X2: ') or line.startswith('# X3: '):
                domain = line.split('[')[1].split(']')[0].split(',')
                bbox.append([float(domain[0]),float(domain[1])])
    return np.array(bbox)

def read_hdf5(hdf5Name,varName,slc=-1,los=0,c=1):
    """
    Loads H5 data from a given file and speficied variable.
    Such variables include:
        'rho' 'prs' 'vx1' 'vx2' 'vx3' 'Bx1' 'Bx2' 'Bx3' 'X' 'Y' 'Z'
    Three-dimensional variables are transposed for alignment with numpy dimensions.
    
    Parameters
    ----------
        hdf5Name : str
            File name (and location) of the .h5 file to load.
        varName : str
            Name of the variable which to be returned.
    
    Returns
    -------
        ret : np.ndarray
           Field values.
    """

    with h5py.File(hdf5Name, 'r') as f:

        if varName in {'rho', 'prs', 'vx1', 'vx2', 'vx3', 'Bx1', 'Bx2', 'Bx3', 'pot'}:
            try:
                ret = np.array(f[list(f.keys())[0]]['vars'][varName])
            except:
                raise KeyError(f'Invalid variable name: {varName}')
        elif varName in ['X', 'Y', 'Z']:
            ret = np.array(f['cell_coords'][varName])
        else:
            raise ValueError(f'Invalid variable name: {varName}')

    if ret.ndim==3:
        ret = ret.T # Transpose if dimension==3
        if slc>-1:
            if los in (0,'x','X'):
                ret = ret[slc][::c,::c]
            elif los in (0,'y','Y'):
                ret = ret[:,slc][::c,::c]
            elif los in (0,'z','Z'):
                ret = ret[:,:,slc][::c,::c]
            else:
                ValueError("Axis must be 0, x, X, 1, y, Y, 2, z, or Z.")
        else:
            ret = ret[::c,::c,::c]
        return ret
    else:
        return ret[::c,::c]

def read_definitions_h(Path:str,param_name:str):
    """
    Retrieve the metadata from the file named definitions.h

    Parameters
    ----------
        Path : str
            Directory which stores all the PLUTO outputs and parameter files.
        param_name : str
            Name of the parameter which to return the value of.
    Returns
    -------
        ret : variable
            If the parameter name is accounted for, it will be formatted appropriately.
            Otherwise, return the data raw.
    """

    with open(os.path.join(Path, "definitions.h"), 'r') as f:
        for line in f.readlines():
            if param_name in line:
                ret = line.split(param_name)[1]
                break
    ValueError(f"The parameter {param_name} was not found in definitions.h.")
    
    if param_name in {'DIMENSIONS'}:
        ret = int(ret)
    elif param_name in {'PHYSICS','EOS','FORCED_TURB'}:
        ret = ret.split(' ')[-1].split('\n')[0]
    elif param_name in {'FORCED_TURB_DECAY', 'FORCED_TURB_ENERGY','FORCED_TURB_KMIN','FORCED_TURB_KMAX'}:
        ret = float(ret)
    elif param_name=='UNIT_DENSITY':
        tmp = ret.split('*')
        ret = 1.
        for _tmp in tmp:
            if 'mp' in _tmp:
                ret *= CGS.mp
            else:
                ret *= float(_tmp)
    elif param_name=='UNIT_LENGTH':
        tmp = ret.split('*')
        ret = 1.
        for _tmp in tmp:
            if 'pc' in _tmp:
                ret *= CGS.pc
            else:
                ret *= float(_tmp)
    elif param_name=='UNIT_VELOCITY':
            ret = float(ret)
    else:
        print(f"Warning: the parameter {param_name} is not formatted for. Returning as is.")

    return ret

def read_pluto_ini(Path,param_name):
    """
    Retrieve the metadata from the file named pluto.ini.

    Parameters
    ----------
        See read_definitions_h
    Returns
    -------
        See read_definitions_h
    """

    with open(os.path.join(Path, "pluto.ini"), 'r') as f:
        for line in f.readlines():
            if param_name in line:
                ret = line.split(param_name)[1]
                break
    ValueError(f"The parameter {param_name} was not found in definitions.h.")

    if param_name in ('X1-grid','X2-grid','X3-grid'):
        ret = [_val for _val in ret.split(' ') if len(_val)>0]

    return ret