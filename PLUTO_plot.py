"""
Script that plots fields from PLUTO.
"""



if __name__=="__main__":

    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser("main")
    parser.add_argument("Path", help="Where the PLUTO outputs are stored.", type=str)
    args = parser.parse_args()
    cluster = args.cluster

    