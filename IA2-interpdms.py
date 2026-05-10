import numpy as np
import h5py, os
from scipy.interpolate import griddata

from utils.enzo import IA2, IA2Data

def smooth_dm(dm,snapshot):
    
    # dm = ia2.get_val('Dark_Matter_Density',snapshot,units='code',c=20)
    
    shape = dm.shape
    xyz = tuple(np.arange(sh) for sh in shape)
    mesh = np.meshgrid(*xyz,indexing='ij')

    where = np.where(dm>0)
    points = tuple(_x[_where] for _x,_where in zip(xyz,where))

    dm_smooth = griddata(points,dm[dm>0],mesh).reshape(shape)
    dm_smooth[np.logical_not(np.isfinite(dm_smooth))] = 0
    dm_smooth = dm_smooth/dm_smooth.sum()*dm.sum()

    with h5py.File(os.path.join(Path,'dm_smooth.hdf5'), 'a') as f:
        f.create_dataset(str(snapshot), data=dm_smooth)

if __name__=="__main__":
    
    Paths = [
        "/media/yange/MyDrive/2024PhDData/EnzoIA2/E14-R/",
        # "/media/yange/MyDrive/2024PhDData/EnzoIA2/E18B-PM/",
        # "/media/yange/MyDrive/2024PhDData/EnzoIA2/E3A-PM/",
        # "/media/yange/MyDrive/2024PhDData/EnzoIA2/E5A-M/",
    ]

    for Path in Paths:
        
        ia2 = IA2Data(Path)

        for s in IA2.snapshots:
            try:
                dm = ia2.get_val('Dark_Matter_Density',s,units='code',c=1)
            except:
                continue
            
            print(f'Interpolating for snapshot={s}.')
            smooth_dm(dm,s)