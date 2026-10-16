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
    get_bbox
        Calculates bounding box of run.
    get_val
        Returns data from specified file with specified variable.
    read_definitions_h
        Return parameter value of speficied parameter.
    read_pluto_ini
        Return parameter value of speficied parameter.

TODO: AMR grid is not accounted for.

Author: Angela Xue
Date: March 2026
Revised: October 2026
"""

import numpy as np
import os, h5py
from utils.units import CGS

class PlutoData:

    def __init__(self,Path):

        self.Path = Path
        self._set_metadata()

    def _set_metadata(self,):
        """
        Retrieve the metadata from the file named definitions.h, pluto.ini and derived quantities.
        """

        # Load all data.XXXX.dbl.h5 names and locations.
        self.hdf5Names = [os.path.join(self.Path,hdf5Name) for hdf5Name in get_hdf5Names(self.Path)]

        self.ndim = read_definitions_h(self.Path,"DIMENSIONS")
        self.PHYSICS = read_definitions_h(self.Path,"PHYSICS")
        self.EOS = read_definitions_h(self.Path,"EOS")
        if self.EOS=="IDEAL":
            self.gamma = 5/3
        try:
            self.FORCED_TURB = read_definitions_h(self.Path,"FORCED_TURB")
        except:
            self.FORCED_TURB = "NO"
        if self.FORCED_TURB=="YES":
            self.FORCED_TURB_DECAY = read_definitions_h(self.Path,"FORCED_TURB_DECAY")
            self.FORCED_TURB_ENERGY = read_definitions_h(self.Path,"FORCED_TURB_ENERGY")
            self.FORCED_TURB_KMIN = read_definitions_h(self.Path,"FORCED_TURB_KMIN")
            self.FORCED_TURB_KMAX = read_definitions_h(self.Path,"FORCED_TURB_KMAX")

        # Units. If not defined, set to default.
        try:
            self.rho_0 = read_definitions_h(self.Path,"UNIT_DENSITY")
        except:
            self.rho_0 = CGS.mp     # proton mass per cm^3
        try:
            self.L_0 = read_definitions_h(self.Path,"UNIT_LENGTH")
        except:
            self.L_0 = CGS.AU       # astronomical unit
        try:
            self.v_0 = read_definitions_h(self.Path,"UNIT_VELOCITY")
        except:
            self.v_0 = 1e5          # 1 km per s
        # Set derived units.
        self.t_0 = self.L_0/self.v_0
        self.p_0 = self.rho_0*self.v_0**2
        self.B_0 = self.v_0*(4*np.pi*self.rho_0)**.5

        # Define shape.
        shape_x, shape_y, shape_z = 1, 1, 1
        shape_x = int(read_pluto_ini(self.Path,"X1-grid")[2])
        shape_y = int(read_pluto_ini(self.Path,"X2-grid")[2])
        shape_z = int(read_pluto_ini(self.Path,"X3-grid")[2])
        self.shape = (shape_x, shape_y, shape_z,)

        self.bbox = get_bbox(self.Path)
        self.varNames = get_varNames(self.Path)

    def print_attr(self,):

        print(f"PHYSICS : {self.PHYSICS}")
        print(f"ndim : {self.ndim}")
        print(f"shape : {self.shape}")
        print(f"bbox : {self.bbox}")
        # print(f"_vars : {self._vars}") #TODO not implemented.
        print(f"EOS : {self.EOS}")
        if self.EOS=="IDEAL":
            print(f"gamma : {self.gamma}")
        print()
        print(f"rho_0 : {self.rho_0} g")
        print(f"L_0 : {self.L_0} cm")
        print(f"v_0 : {self.v_0} cm/s")
        print(f"t_0 : {self.t_0} s")
        print(f"p_0 : {self.p_0} cgs")
        print(f"B_0 : {self.B_0} cgs")
        print()
        print(f"FORCED_TURB : {self.FORCED_TURB}")
        if self.FORCED_TURB=="YES":
            print(f"FORCED_TURB_DECAY : {self.FORCED_TURB_DECAY}")
            print(f"FORCED_TURB_ENERGY : {self.FORCED_TURB_ENERGY}")

    def get_val(self,varName,snapshot,units="code",**kwargs):
            """
            Load and return field values of a single field from a specified snapshot.
            See get_val for details.
            """
    
            _const = 1.
            if units.lower()=="cgs":
                if varName=="rho":
                    _const *= self.rho_0
                elif "vx" in varName:
                    _const *= self.v_0
                elif "Bx" in varName:
                    _const *= self.B_0
            elif units.lower()=="code":
                _const = 1.
            else:
                raise ValueError(f"Units {units} is not supported.")
    
            ret = get_val(self.Path,varName,snapshot,**kwargs) * _const
            return ret

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
    
    def get_times(self,snapshots="all",units="code"):
        """ Return time in either code or CGS units. """

        if units=='code':
            _const = 1.
        elif units=='cgs':
            _const = self.t_0
        else:
            print(f"The unit {units} is not supported; defaulting to code units.")
            _const = 1.

        times = _const*get_times(self.Path)
        if type(snapshots)==int:
            times = times[snapshots]
        elif type(snapshots)==list:
            times = np.array([times[i] for i in snapshots])
        elif snapshots=="all":
            times = times
        else:
            ValueError(f"snapshots = {snapshots} is not supported.")
        return times

    def get_snapshots(self,):
        """ Read file names in the path to retrieve list of snapshots. """
        fNames = []
        for fName in os.listdir(self.Path):
            if fName.startswith("data") and fName.endswith("dbl.h5"):
                fNames.append(fName)
        snapshots = np.array([int(fName.split(".")[1]) for fName in fNames])
        idx = np.argsort(snapshots)
        snapshots = snapshots[idx]
        return snapshots        

    def get_dxdydz(self,units="code"):
        """ Get cell size. """

        if units=='code':
            _const = 1.
        elif units=='cgs':
            _const = self.L_0
        else:
            print(f"The unit {units} is not supported; defaulting to code units.")
            _const = 1.

        dxdydz = np.array([np.diff(_b)/_s for _b,_s in zip(self.bbox,self.shape)])
        dxdydz *= _const
        return dxdydz


def get_hdf5Names(Path:str):
    """
    Load all HDF5 file names in a directory in ascending order.
    A maximum of 1000 output files are produced by PLUTO.
    """
    f_names = [f_name for f_name in os.listdir(Path) if f_name.endswith(".dbl.h5")]
    f_nums = np.sort([int(f_name[5:9]) for f_name in f_names])
    sorted_hdf5Names = [f"data.{("000"+str(f_num))[-4:]}.dbl.h5" for f_num in f_nums]
    return sorted_hdf5Names

def get_varNames(Path:str):
    """ Retrieve field variable names. """
    fName = get_hdf5Names(Path)[0]
    with h5py.File(os.path.join(Path,fName)) as f:
        varNames = list(f["Timestep_0"]["vars"].keys())
    return varNames

def get_times(Path:str):
    """ Get all times from the file "dbl.h5.out". Code units. """
    return np.loadtxt(os.path.join(Path,"dbl.h5.out"), usecols=1)

def get_bbox(Path:str):
    """ Retrive bounding box. Code units. """
    bbox = []
    with open(os.path.join(Path,"grid.out")) as f:
        for line in f.readlines()[:20]: # only need to parse the top of the file.
            if line.startswith("# X1: ") or line.startswith("# X2: ") or line.startswith("# X3: "):
                domain = line.split("[")[1].split("]")[0].split(",")
                bbox.append([float(domain[0]),float(domain[1])])
    return np.array(bbox)

def get_val(Path,varName,snapshot,c=1):
    """
    Loads H5 data from a given file and speficied variable.
    Such variables include:
        'rho' 'prs' 'vx1' 'vx2' 'vx3' 'Bx1' 'Bx2' 'Bx3' 'X' 'Y' 'Z'
    Three-dimensional variables are transposed for alignment with numpy dimensions.
    
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

    snapshot_str = ("0000"+str(snapshot))[-4:]
    hdf5Name = os.path.join(Path,"data."+snapshot_str+".dbl.h5")
    with h5py.File(hdf5Name, "r") as f:
        if varName in {"rho","prs","vx1","vx2","vx3","Bx1","Bx2","Bx3","pot"}:
            try:
                ret = np.array(f[list(f.keys())[0]]["vars"][varName])
            except:
                raise ValueError(f"Variable {varName} is not recognised.")
        elif varName in ["X", "Y", "Z"]:
            ret = np.array(f["cell_coords"][varName])
        else:
            raise ValueError(f"Variable {varName} is not recognised.")

    if ret.ndim==3:
        ret = ret.T # Transpose if dimension==3
        return ret[::c,::c,::c]
    else:
        return ret[::c,::c]

def read_definitions_h(Path,paramName):
    """
    Retrieve the metadata from the file named definitions.h

    Parameters
    ----------
        Path : str
            Directory which stores all the PLUTO outputs and parameter files.
        paramName : str
            Name of the parameter which to return the value of.
    Returns
    -------
        ret : variable
            If the parameter name is accounted for, it will be formatted appropriately.
            Otherwise, return the data raw.
    """

    ret = None
    with open(os.path.join(Path, "definitions.h"), "r") as f:
        for line in f.readlines():
            line = line.removesuffix("\n")
            if paramName in line:
                ret = line.split(paramName)[1]
                break
    try:
        assert type(ret)==str
    except:
        raise ValueError(f"The parameter {paramName} was not found in definitions.h.")
    
    if paramName in {"DIMENSIONS"}:
        ret = int(ret)
    elif paramName in {"PHYSICS","EOS","FORCED_TURB"}:
        ret = ret.split(" ")[-1].split("\n")[0]
    elif paramName in {"FORCED_TURB_DECAY", "FORCED_TURB_ENERGY","FORCED_TURB_KMIN","FORCED_TURB_KMAX"}:
        ret = float(ret)
    elif paramName=="UNIT_DENSITY":
        tmp = ret.split("*")
        ret = 1.
        for _tmp in tmp:
            if "mp" in _tmp:
                ret *= CGS.mp
            else:
                ret *= float(_tmp)
    elif paramName=="UNIT_LENGTH":
        tmp = ret.split("*")
        ret = 1.
        for _tmp in tmp:
            if "pc" in _tmp:
                ret *= CGS.pc
            else:
                ret *= float(_tmp)
    elif paramName=="UNIT_VELOCITY":
            ret = float(ret)
    else:
        print(f"Warning: the parameter {paramName} is not formatted for. Returning as is.")

    return ret

def read_pluto_ini(Path,paramName):
    """
    Retrieve the metadata from the file named pluto.ini.
    Special considerations are given for reading the grid.

    Parameters
    ----------
        See read_definitions_h.
    Returns
    -------
        See read_definitions_h.
    """

    ret = None
    with open(os.path.join(Path, "pluto.ini"), "r") as f:
        for line in f.readlines():
            line = line.removesuffix("\n")
            if line.startswith(paramName):
                ret = line.split(paramName)[1]
                break
    try:
        assert type(ret)==str
    except:
        raise ValueError(f"The parameter {paramName} was not found in definitions.h.")
    
    ret = [_val for _val in ret.split(" ") if len(_val)>0]

    return ret