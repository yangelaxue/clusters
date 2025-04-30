"""
Functions that help with reading and loading data from PLUTO simulations.

Author: Angela Xue
"""

#%% Imports

import numpy as np
import os, h5py


#%% Functions that help with loading outout data.

def get_times(output_dir):
    """
    Get the times from the file 'dbl.h5.out'.   
    """

    return np.loadtxt(os.path.join(output_dir,'dbl.h5.out'), usecols=1)

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

def load_HDF5data(f_name, *var):
    """
    Loads H5 data from a given file and speficied variable.
    Such variables include:
        'rho' 'prs' 'vx1' 'vx2' 'vx3' 'Bx1' 'Bx2' 'Bx3' 
    
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

    return out

def get_units(f_name):
    """
    Retrieve the scaling units from the file named definitions.h
    These scaling units are L_0, v_0 and rho_0, all listed in cgs units.
    If these units aren't found, default to 1.

    Parameters
    ----------
        f_name : str
            File name (and location) of definitions.h file to search within.
    
    Returns
    -------
        L_0, v_0, rho_0 : float
            The units scaling units used by PLUTO.
    """

    L_bool, v_bool, rho_bool = 0, 0, 0

    with open(f_name, 'r') as f:
        for line in f.readlines():
            if line.startswith('#define UNIT_LENGTH'):
                L_0 = float(line.split('//')[0].split('UNIT_LENGTH')[1])
                L_bool = 1
            elif line.startswith('#define UNIT_VELOCITY'):
                v_0 = float(line.split('//')[0].split('UNIT_VELOCITY')[1])
                v_bool = 1
            elif line.startswith('#define UNIT_DENSITY'):
                rho_0 = float(line.split('//')[0].split('UNIT_DENSITY')[1])
                rho_bool = 1

    L_0 = L_0 if L_bool else 1.
    v_0 = v_0 if v_bool else 1.
    rho_0 = rho_0 if rho_bool else 1.

    return L_0, v_0, rho_0