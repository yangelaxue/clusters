"""
Read data from different simulations into the same format.

Programs supported:
    - PLUTO
"""

import h5py
import os
import numpy as np

#%% PLUTO

def get_H5_filenames(f_dir):
    """
    Returns all of the data file names in order.
    """

    f_idxes = np.sort([int(_fname.split(".")[1]) for _fname in os.listdir(f_dir) if _fname.endswith(".dbl.h5")])
    f_names = [f"data.{('000'+str(_f_idx))[-4:]}.dbl.h5" for _f_idx in f_idxes]
    
    return f_names

def read_H5(f_dir, *keys):
    """
    Reads and returns values from PLUTO H5 files. This includes field values as well as coordinate values.
    
    Currently only supports static grid outputs.
    
    ----------
    Parameters
    ----------
    f_dir : str     Directory and file name of the file to be read.
    *keys : str     Keys intended of the values to be returned.

    -------
    Returns
    -------
    vals : list     list of values to be 
    """

    # Checks before proceeding to main function.
    for key in keys:
        assert type(key)==str, "Keys must be type str"

    # START
    file = h5py.File(f_dir, 'r')

    vals = []

    for key in keys:
        if key in ['X', 'Y', 'Z']:
            vals.append(file['cell_coords'][key][...])
        else:
            vals.append(file[list(file.keys())[0]]['vars'][key][...])

    return vals
    
    



    




