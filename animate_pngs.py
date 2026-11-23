"""
Script that animates a list of PNGs.

Author: Angela Xue
Date: Nov 2026
"""

if __name__ == "__main__":

    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Path with PNGs to animate.", type=str)
    parser.add_argument("fName", help="Name prefix of PNGs to animate.", type=str)
    args = parser.parse_args()
    Path = args.Path
    fName = args.fName

    import numpy as np
    import os
    from utils.visualise import animate_pngs

    fNames = np.array([f for f in os.listdir(Path) if f.startswith(fName)])
    idxes = [int(f.split(".")[0].split("_")[-1]) for f in fNames]
    argsort = np.argsort(idxes)
    print(argsort)
    fNames = fNames[argsort]

    fNames = [os.path.join(Path, f) for f in fNames]
    animate_pngs(fNames,figsize=(10,10),interval=300,)

    