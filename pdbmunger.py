#!/usr/bin/env python
# Copyright (C) 2021 Daniel Asarnow
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import argparse
import pandas as pd
import numpy as np
import sys
import gzip
from biopandas.pdb import PandasPdb


def main(args):
    pdb = PandasPdb().read_pdb(args.input)
    if args.renumber_atoms:
        pdb.df["ATOM"]["atom_number"] = pdb.df["ATOM"].index.values + 1
        pdb.df["HETATM"]["atom_number"] = pdb.df["HETATM"].index.values + pdb.df["ATOM"].shape[0]
    if args.occupancy is not None:
        pdb.df["ATOM"].loc[:, "occupancy"] = np.float(args.occupancy)
        if "HETATM" in pdb.df:
            pdb.df["HETATM"].loc[:, "occupancy"] = np.float(args.occupancy)
    if args.bfactor is not None:
        bfact = PandasPdb().read_pdb(args.bfactor)
        bfactatoms = bfact.df["ATOM"].set_index(["chain_id", "residue_number", "atom_name"])
        atoms = pdb.df["ATOM"].set_index(["chain_id", "residue_number", "atom_name"])
        print("Updating %d protein atoms" % bfactatoms.loc[atoms.index, "b_factor"].shape[0])
        atoms.loc[:, "b_factor"] = bfactatoms.loc[atoms.index, "b_factor"]
        pdb.df["ATOM"] = atoms.reset_index()[pdb.df["ATOM"].columns]
        if "HETATM" in pdb.df:
            hetatoms = pdb.df["HETATM"].set_index(["chain_id", "residue_number", "atom_name"])
            print("Updating %d heteroatoms" % bfactatoms.loc[hetatoms.index, "b_factor"].shape[0])
            hetatoms.loc[:, "b_factor"] = bfactatoms.loc[hetatoms.index, "b_factor"]
            pdb.df["HETATM"] = hetatoms.reset_index()[pdb.df["HETATM"].columns]
    pdb.to_pdb(args.output)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--renumber-atoms", "-ra", action="store_true",
            help="Renumber ATOM and HETATM records from 1 (WILL DELETE CONECT RECORDS)")
    parser.add_argument("--bfactor", "-b", help="Update b-factors from source PDB or to single value")
    parser.add_argument("--occupancy", "-o", help="Update occupancies from source PDB or to single value")
    sys.exit(main(parser.parse_args()))

