#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import sys
import gzip
from biopandas.pdb import PandasPdb


def main(args):
    pdb = PandasPdb().read_pdb(args.input)
    if args.bfactor is not None:
        bfact = PandasPdb().read_pdb(args.bfactor)
        bfactatoms = bfact.df["ATOM"].set_index(["chain_id", "residue_number", "atom_name"])
        atoms = pdb.df["ATOM"].set_index(["chain_id", "residue_number", "atom_name"])
        atoms.loc[:, "b_factor"] = bfactidx.loc[atoms.index, "b_factor"]
        hetatoms = pdb.df["HETATM"].set_index(["chain_id", "residue_number", "atom_name"])
        hetatoms.loc[:, "b_factor"] = bfactatoms.loc[hetatoms.index, "b_factor"]
        pdb.df["ATOM"] = atoms
        pdb.df["HETATOM"] = hetatoms
    pdb.to_pdb(args.output)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--bfactor", "-b", help="Update b-factors from source PDB or to single value")
    sys.exit(main(parser.parse_args()))

