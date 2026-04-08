"""
Calculate and save cluster centres from the IA2 clusters.
"""

if __name__=="__main__":

    from utils.enzo import IA2Data, IA2

    Paths = [
        "/nesi/nobackup/vuw04655/IA2/E14-R",
        "/nesi/nobackup/vuw04655/IA2/E18B-PM",
        "/nesi/nobackup/vuw04655/IA2/E3A-PM",
        "/nesi/nobackup/vuw04655/IA2/E5A-M",
    ]

    for Path in Paths:

        ia2 = IA2Data(Path)

        for s in IA2.snapshots:
            try:
                ia2.get_centre(s,_c=5)
            except:
                continue

