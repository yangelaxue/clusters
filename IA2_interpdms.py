import numpy as np
import h5py, os
from scipy.interpolate import griddata

def smooth_dm(dm):
    
    shape = dm.shape
    xyz = tuple(np.arange(sh) for sh in shape)
    mesh = np.meshgrid(*xyz,indexing='ij')

    where = np.where(dm>0)
    points = tuple(_x[_where] for _x,_where in zip(xyz,where))

    dm_smooth = griddata(points,dm[dm>0],mesh).reshape(shape)
    dm_smooth[np.logical_not(np.isfinite(dm_smooth))] = 0
    dm_smooth = dm_smooth/dm_smooth.sum()*dm.sum()

    return dm_smooth

def loop_success(sns):

    dms = ia2.gen_vals('Dark_Matter_Density',snapshots=sns,units='code',c=c)

    idx = 0
    for s in sns:
        
        if os.path.exists(saveName):
            with h5py.File(saveName, 'r') as f:
                if str(s) in f.keys():
                    print(f"Interpolated DM for snapshot {s} is already calculated.")
                    idx += 1
                    continue
        try:
            dm = next(dms)
        except:
            print(f"Error occured for snapshot={s}. Skipping.")
            if len(sns)>idx+1:
                return sns[idx+1:]
            else:
                print("Success.")
                return []
        dm_smooth = smooth_dm(dm)
        print(f'Saving interpolation for snapshot={s}.')
        with h5py.File(saveName, 'a') as f:
            f.create_dataset(str(s), data=dm_smooth)
        idx += 1
    return []

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
    snapshots = args.snapshots if args.snapshots else 'all'
    if snapshots=='all':
        snapshots = list(IA2.snapshots)
    else:
        snapshots = args.snapshots.removeprefix('[').removesuffix(']')
        snapshots = [int(s) for s in snapshots.split(',')]
    

    saveName = os.path.join(Path,f'dm_smooth_{c}.hdf5')
    ia2 = IA2Data(Path)
    
    # If snapshot has already been interpolated, remove it from the snapshots.
    if os.path.exists(saveName):
        tmp = []
        with h5py.File(saveName, 'r') as f:
            for s in snapshots:
                if str(s) not in f.keys():
                    tmp.append(s)
                else:
                    print(f"Already interpolated snapshot {s}.")
        snapshots = tmp
    # If corresponding file does not exist, remove snapshot from snapshots.
    tmp = []
    fNames = ia2.get_fnames('dt')['dt']
    for s in snapshots:
        for fName in fNames:
            if fName.endswith(('000'+str(s))[-3:]):
                tmp.append(s)
    snapshots = tmp

    sns = snapshots
    while sns:
        sns = loop_success(sns)
    