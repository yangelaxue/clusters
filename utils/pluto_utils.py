"""
Functions that help with reading and loading data from PLUTO simulations.

Author: Angela Xue
"""

#%% Imports

import numpy as np
import os, h5py

#%% Create class

class PlutoUnits:
    """
    PLUTO Code uses cgs units by default.
    Helps with:
        - calculating other unit quantities from the given cgs units,
        - Converting from cgs units to other units.
    """

    # List important units in cgs units.

    CONST_AU = 1.49597892e13        # Astronomical unit
    CONST_mp = 1.67262171e-24       # Proton mass
    CONST_c = 2.99792458e10         # Light speed
    CONST_e = 4.80320425e-10        # Elementary proton charge
    CONST_eV = 1.602176463158e-12   # Electron Volt in erg.
    CONST_pc = 3.0856775807e18      # Parcec in cm

    def set_units_from_file(self,f_dir):
        """
        Retrieve the scaling units from the file named definitions.h
        These scaling units are rho_0, L_0 and v_0, all listed in cgs units.
        If these units aren't found, use the default settings used by PLUTOm which are:
            rho_0   = 1 proton mass per cm^3    = 1.67262171e-24    g/cm^3
            L_0     = 1 AU                      = 1.49597892e13     cm
            v_0     = 1 km/s                    = 100,000           cm/s

        Parameters
        ----------
            f_name : str
                Location of definitions.h file to search within.
        """

        # Append difinitions.h file name to directory if need be.
        if not f_dir.endswith("definitions.h"):
            f_dir = os.path.join(f_dir, "definitions.h")

        # Set units to default PLUTO values.
        rho_0 = self.CONST_mp
        L_0 = self.CONST_AU
        v_0 = 1e5

        # Loop over all lines, redefine units if found.
        with open(f_dir, 'r') as f:
            for line in f.readlines():
                if 'UNIT_DENSITY' in line:
                    rho_0 = line.split('UNIT_DENSITY')[1]
                    if 'CONST_mp' in rho_0:
                        rho_0 = self.CONST_mp
                    else:
                        rho_0 = float(rho_0)
                elif 'UNIT_LENGTH' in line:
                    L_0 = float(line.split('UNIT_LENGTH')[1])
                elif 'UNIT_VELOCITY' in line:
                    v_0 = float(line.split('UNIT_VELOCITY')[1])

        self.rho_0, self.L_0, self.v_0 = rho_0, L_0, v_0

    def set_T_0(self,):
        self.T_0 = self.L_0/self.v_0

class PlutoData:

    def __init__(self, output_dir:str):

        self.output_dir = output_dir
        # load_metadata()
        self.load_metadata()

        # Set units.
        self.Units = PlutoUnits()
        self.Units.set_units_from_file(output_dir)
        self.Units.set_T_0()

    def load_metadata(self):
        """
        Retrieve the metadata from the file named definitions.h
        """

        # Loop over all lines, redefine units if found.
        with open(os.path.join(self.output_dir, "definitions.h"), 'r') as f:
            for line in f.readlines():
                if 'DIMENSIONS' in line:
                    self.ndim = int(line.split('DIMENSIONS')[1])
    
    def get_HDF5names(self,):
        """
        Return all data.XXXX.dbl.h5 file names.
        """
        return get_HDF5names(self.output_dir)
    
    def load_times(self,units:str):
        """
        Return all time slices, in code units.
        """

        assert units in ['code','physical'], f"Units must be 'code' or 'physical', not {units}."
        if units=='code':
            return load_times(self.output_dir)
        else:
            return load_times(self.output_dir) * self.Units.T_0
    
    def load_domain_from_HDF5(self,units:str):
        """
        Return domain in physical units.
        """

        assert units in ['code','physical'], f"Units must be 'code' or 'physical', not {units}."
        var = ['X']
        if self.ndim>1:
            var.append('Y')
        if self.ndim>2:
            var.append('Z')
        ret = load_data_from_HDF5(os.path.join(self.output_dir,self.get_HDF5names()[0]), *var)

        if units=='physical':
            for i, _ret in enumerate(ret):
                ret[i] = _ret  * self._get_conversion_from_var('X')

        return ret
    
    def load_data_from_HDF5(self, var:str, units:str, timeslices='all'):
        """
        Return variables in physical units.

        Parameters
        ----------
            var : str
                Name of the variable to be loaded.
            units: str
                Whether to return the variable in 'code' or 'physics' units.
            timeslices : str, arraylike
                List of timeslices to load variables from. If 'all', return all timeslices.

        Returns
        
        """

        assert units in ['code','physical'], f"Units must be 'code' or 'physical', not {units}."
        
        timeslices = [timeslices] if type(timeslices)==int else timeslices

        f_names = self.get_HDF5names()
        if timeslices!='all':
            assert type(timeslices)==list, "timeslices needs to be a list of integers"
            f_names = [self.get_HDF5names()[timeslice] for timeslice in timeslices]

        if units=='code':
            return [load_data_from_HDF5(os.path.join(self.output_dir,f_name), var) for f_name in f_names]
        else:
            return [load_data_from_HDF5(os.path.join(self.output_dir,f_name), var)*self._get_conversion_from_var(var) for f_name in f_names]
            
    def _get_conversion_from_var(self,var):
        """
        Given some variable name, return the conversion from code units to cgs units.
        Used for converting units of HDF5 data.
        """
        if var=='rho':
            return self.Units.rho_0
        elif var=='prs':
            return self.Units.rho_0*self.Units.v_0**2
        elif var in ['vx1','vx2','vx3']:
            return self.Units.v_0
        elif var in ['Bx1','Bx2','Bx3']:
            return (4*np.pi*self.Units.rho_0*self.Units.v_0**2)**.5
        elif var in ['X','Y','Z']:
            return self.Units.L_0
        else:
            print(f"Variable {var} not recognised - defaulting to 1.0...\n")
            return 1.0


#%% Functions that help with loading outout data.

def get_HDF5names(output_dir):
    """
    Load all HDF5 file names in a directory in ascending order.
    A maximum of 1000 output files are produced by PLUTO.

    Parameters
    ----------
        output_dir : str
            Directory which to list all HDF5 files.
    
    Returns
    -------
        sorted_f_names : str
            A list of all HDF5 file names in the order of their names..
    """

    f_names = [f_name for f_name in os.listdir(output_dir) if f_name.endswith('.dbl.h5')]
    f_nums = np.sort([int(f_name[5:9]) for f_name in f_names])
    sorted_f_names = [f'data.{('000'+str(f_num))[-4:]}.dbl.h5' for f_num in f_nums]

    return sorted_f_names

def load_times(output_dir):
    """
    Get the times from the file 'dbl.h5.out'.   
    """

    return np.loadtxt(os.path.join(output_dir,'dbl.h5.out'), usecols=1)

def load_data_from_HDF5(f_name, *var):
    """
    Loads H5 data from a given file and speficied variable.
    Such variables include:
        'rho' 'prs' 'vx1' 'vx2' 'vx3' 'Bx1' 'Bx2' 'Bx3' 'X' 'Y' 'Z'
    
    Parameters
    ----------
        f_name : str
            File name (and location) of the .h5 file to load.
        *var : str
            Names of the variables which to be returned.
    
    Returns
    -------
        out : list
            List of the grid of variables returned in the same order as var.
    """

    out = []

    with h5py.File(f_name, 'r') as f:

        for _var in var:
            if _var in ['rho', 'prs', 'vx1', 'vx2', 'vx3', 'Bx1', 'Bx2', 'Bx3']:
                out.append(np.array(f[list(f.keys())[0]]['vars'][_var]))
            elif _var in ['X', 'Y', 'Z']:
                out.append(np.array(f['cell_coords'][_var]))
            else:
                print(f'Invalid variable name: {_var}')
                out.append(0) # Return 0 to ensure all variables are in the same order.

    out = out if len(out)>1 else out[0]
    return out