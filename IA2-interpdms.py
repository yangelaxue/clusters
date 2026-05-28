import numpy as np
import h5py, os
from scipy.interpolate import griddata

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

    import argparse, os
    from utils.enzo import IA2, IA2Data

    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path to PLUTO data.", type=str)
    parser.add_argument("-c", help="How small a subset to plot.", type=int,required=True)
    parser.add_argument("-snapshots", help="Which snapshots to animate.",required=False)

    args = parser.parse_args()
    assert os.path.exists(args.Path), "First argument MUST be an existing path."
    Path = args.Path
    c = args.c
    if not args.snapshots:
        snapshots = IA2.snapshots
    elif args.snapshots!='all':
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]
    
    ia2 = IA2Data(Path)
    tmp = []
    fNames = ia2.get_fnames('b')
    for s in snapshots:
        for fName in fNames:
            if fName.endswith(str(s)):
                tmp.append(s)
                continue
    snapshots = tmp

    dms = ia2.gen_vals('Dark_Matter_Density',snapshots=snapshots,units='code',c=c)
    for s,dm in zip(snapshots,dms):
        print(f'Interpolating for snapshot={s}.')
        smooth_dm(dm,s)
