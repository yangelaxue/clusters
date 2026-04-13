"""
Calculate and save cluster centres from the IA2 clusters.
"""

if __name__=="__main__":

    import numpy as np
    import os
    from utils.enzo import IA2Data, IA2

    Paths = [
        "/nesi/nobackup/vuw04655/IA2/E14-R",
        "/nesi/nobackup/vuw04655/IA2/E18B-PM",
        "/nesi/nobackup/vuw04655/IA2/E3A-PM",
        "/nesi/nobackup/vuw04655/IA2/E5A-M",
    ]

    c = 1

    for Path in Paths:

        ia2 = IA2Data(Path)
        rho_cs = []

        for s in IA2.snapshots:
            try:
                centre = ia2.get_centre(s,_c=c)
            except:
                continue
            
            rho_cs.append([s])
            rho_c = ia2.get_val('Density',s,c=c,slc=int(centre[0]))[int(centre[1]),int(centre[2])]
            rho_cs[-1].append(rho_c)
        
        np.savetxt(os.path.join(Path,'rho_c.txt'),rho_cs,header='snapshot rho_c')
            

