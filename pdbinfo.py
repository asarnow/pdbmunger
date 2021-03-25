#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import sys
import gzip
from biopandas.pdb import PandasPdb


def main(args):
    pdb = PandasPdb().read_pdb(args.input)
    print("Protein atoms: %d (%d non-hydrogen)" %
            (pdb.df["ATOM"].shape[0], pdb.get('heavy').shape[0]))
    print("Hetero atoms: %d (%d non-hydrogen)" %
            (pdb.df["HETATM"].shape[0], pdb.df["HETATM"][(pdb.df["HETATM"]['element_symbol'] != 'H')].shape[0]))
    print("Coordinate centroid: %f, %f, %f" %
            (pdb.df["ATOM"][["x_coord"]].mean(), pdb.df["ATOM"][["y_coord"]].mean(), pdb.df["ATOM"][["z_coord"]].mean()))
    print("B-factors (min/max/mean)\n\tProtein: %f / %f / %f\n\tHetero: %f / %f / %f" %
            (pdb.df["ATOM"]["b_factor"].min(), pdb.df["ATOM"]["b_factor"].max(), pdb.df["ATOM"]["b_factor"].mean(),
             pdb.df["HETATM"]["b_factor"].min(), pdb.df["HETATM"]["b_factor"].max(), pdb.df["HETATM"]["b_factor"].mean()))
    print("Occupancies (min/max/mean)\n\tProtein: %f / %f / %f\n\tHetero: %f / %f / %f" % 
            (pdb.df["ATOM"]["occupancy"].min(), pdb.df["ATOM"]["occupancy"].max(), pdb.df["ATOM"]["occupancy"].mean(),
            pdb.df["HETATM"]["occupancy"].min(), pdb.df["HETATM"]["occupancy"].max(), pdb.df["HETATM"]["occupancy"].mean()))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    sys.exit(main(parser.parse_args()))

