"""
Functions used to manipulate volumes of field values

Author: Angela Xue
"""

import numpy as np

#%%

def shift_volume(arr, idxes):
    """
    Redefines a volume of data arr to center a chosen point idxes. It does this by
        slicing the arr along each axes determined by the desired center and glueing
        the opposite faces together, i.e. making use of the periodic spacial boundary conditions.

    Parameters
    ----------
    arr : np.ndarray
        Array of data which is to be shifted and glued.
    idxes : tuple
        Indices of the point which is to be moved to the center of the volume.

    Returns
    -------
    arr_sh : np.ndarray
        Shifted values.
    """
    
    shape = arr.shape
    c_idxes = tuple(shp//2 for shp in shape)

    # Set up shifted values
    arr_sh = arr.copy()

    for i, (c_idx, idx) in enumerate(zip(c_idxes, idxes)):
        
        # Initiate slices.
        slices_l = [slice(0,shape_x) for shape_x in shape]
        slices_r = [slice(0,shape_x) for shape_x in shape]

        # Set up slices to cut and glue to.
        shift = (idx-c_idx)%shape[i]
        slices_l[i] = slice(0,shift)
        slices_r[i] = slice(shift,None)
        
        # Shift array.
        _arr_sh = arr_sh.copy()
        arr_sh = np.concatenate((_arr_sh,_arr_sh[tuple(slices_l)]),axis=i)
        arr_sh = arr_sh[tuple(slices_r)]
            
    return arr_sh