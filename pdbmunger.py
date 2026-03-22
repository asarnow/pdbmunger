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
import sys
import click
import pandas as pd
from biopandas.pdb import PandasPdb


AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    "X": "UNK",
}

# Record types that belong after ATOM/HETATM in canonical PDB order.
# TER is excluded: it is regenerated between chain blocks, not treated as a trailer.
_TRAILER_RECORDS = frozenset({"ANISOU", "ENDMDL", "END", "CONECT", "MASTER"})


def _sort_residues(df, chain_ids=None):
    """Sort residue entries by residue number within each chain, preserving original
    chain order and atom order within each residue.
    If chain_ids is given, only sort those chains; others are left in place."""
    chain_order = df["chain_id"].unique()
    parts = []
    for cid in chain_order:
        sub = df[df["chain_id"] == cid].copy()
        if chain_ids is None or cid in chain_ids:
            sub = sub.sort_values(["residue_number", "atom_number"])
        parts.append(sub)
    return pd.concat(parts).reset_index(drop=True)


def _renumber_residues(df, start, chain_ids=None):
    """Renumber residues in df per chain sequentially from start.
    If chain_ids is given, only renumber those chains.
    Insertion codes are cleared for modified chains."""
    df = df.copy()
    for chain_id in df["chain_id"].unique():
        if chain_ids is not None and chain_id not in chain_ids:
            continue
        mask = df["chain_id"] == chain_id
        old_nums = list(dict.fromkeys(df.loc[mask, "residue_number"]))
        mapping = {old: new for new, old in enumerate(old_nums, start=start)}
        df.loc[mask, "residue_number"] = df.loc[mask, "residue_number"].map(mapping)
        df.loc[mask, "insertion"] = ""
    return df


def _get_target_chains(pdb, chain_opt):
    """Return set of ATOM chain IDs to operate on, or None for all chains.

    Uses amino3to1 to build per-chain sequences for identity comparison.
    - No chain_opt and all ATOM chains identical: return all chain IDs.
    - No chain_opt and multiple distinct groups: return None (apply to all).
    - chain_opt given: expand each specified chain to all ATOM chains with the
      same sequence; include specified chains not found in ATOM (HETATM-only).
    """
    seq_df = pdb.amino3to1(record="ATOM", residue_col="residue_name", fillna="?")
    sequences = seq_df.groupby("chain_id", sort=False)["residue_name"].apply("".join).to_dict()
    if chain_opt is None:
        if len(set(sequences.values())) == 1:
            return set(sequences.keys())
        return None
    specified = {c.strip() for c in chain_opt.split(",")}
    target = set()
    for cid in specified:
        if cid in sequences:
            seq = sequences[cid]
            target.update(c for c, s in sequences.items() if s == seq)
        else:
            target.add(cid)  # may be HETATM-only; handled separately
    return target


def _resolve_hetatm_collisions(hetatm_df, atom_df):
    """Ensure HETATM residue numbers fall outside the ATOM residue range in mixed chains.

    For each chain present in both ATOM and HETATM, if any HETATM residue number
    falls within [min_atom_res, max_atom_res], renumber all HETATM residues for that
    chain sequentially from max_atom_res + 1 (preserving their original order).
    """
    hetatm_df = hetatm_df.copy()
    mixed_chains = set(hetatm_df["chain_id"]) & set(atom_df["chain_id"])
    for chain_id in mixed_chains:
        atom_mask = atom_df["chain_id"] == chain_id
        hetatm_mask = hetatm_df["chain_id"] == chain_id
        atom_min = atom_df.loc[atom_mask, "residue_number"].min()
        atom_max = atom_df.loc[atom_mask, "residue_number"].max()
        hetatm_res = hetatm_df.loc[hetatm_mask, "residue_number"]
        if ((hetatm_res >= atom_min) & (hetatm_res <= atom_max)).any():
            old_nums = list(dict.fromkeys(hetatm_df.loc[hetatm_mask, "residue_number"]))
            mapping = {old: new for new, old in enumerate(old_nums, start=int(atom_max) + 1)}
            hetatm_df.loc[hetatm_mask, "residue_number"] = (
                hetatm_df.loc[hetatm_mask, "residue_number"].map(mapping))
    return hetatm_df


def _parse_fasta(path):
    """Parse a FASTA file and return {record_name: sequence} in file order."""
    chains = {}
    current_id = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    chains[current_id] = "".join(seq_parts)
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if current_id is not None:
        chains[current_id] = "".join(seq_parts)
    return chains


def _seqres_entries(chains):
    """Return list of SEQRES entry strings (columns 7-80) for all chains.

    Each chain maps a chain ID to a single-letter amino acid sequence.
    13 residues are written per record line.
    """
    entries = []
    for chain_id, sequence in chains.items():
        residues = []
        for c in sequence:
            if c not in AA3:
                click.echo("Warning: unknown residue '%s', substituting UNK" % c, err=True)
                residues.append("UNK")
            else:
                residues.append(AA3[c])
        num_res = len(residues)
        chunks = [residues[i:i + 13] for i in range(0, len(residues), 13)]
        for ser_num, chunk in enumerate(chunks, 1):
            res_str = " ".join("%3s" % r for r in chunk)
            entries.append(" %3d %s %4d  %s" % (ser_num, chain_id, num_res, res_str))
    return entries


def _ter_entry(last_atom):
    """Format the entry string (columns 7-80) for a TER record terminating a chain."""
    icode = str(last_atom["insertion"]).strip()
    return "%5d      %3s %s%4d%s" % (
        int(last_atom["atom_number"]) + 1,
        last_atom["residue_name"].strip(),
        last_atom["chain_id"],
        int(last_atom["residue_number"]),
        icode,
    )


def _rebuild_line_idx(pdb, new_seqres=None):
    """Assign line_idx to all records in canonical PDB order:
    header OTHERS -> SEQRES -> per-chain (ATOM + TER) -> HETATM -> trailer OTHERS.

    Existing TER records are dropped and regenerated from current ATOM data.
    Within each group the original relative order is preserved.
    If new_seqres is provided, existing SEQRES records are replaced with it.
    """
    others = pdb.df["OTHERS"]
    names = others["record_name"].str.strip()
    is_seqres = names == "SEQRES"
    is_ter = names == "TER"
    is_trailer = names.isin(_TRAILER_RECORDS)

    header = others[~is_seqres & ~is_ter & ~is_trailer].sort_values("line_idx").copy()
    trailer = others[is_trailer].sort_values("line_idx").copy()

    idx = 0
    header["line_idx"] = range(idx, idx + len(header))
    idx += len(header)

    if new_seqres is not None:
        seqres_df = pd.DataFrame({
            "record_name": ["SEQRES"] * len(new_seqres),
            "entry": new_seqres,
            "line_idx": range(idx, idx + len(new_seqres)),
        })
    else:
        seqres_df = others[is_seqres].sort_values("line_idx").copy()
        seqres_df["line_idx"] = range(idx, idx + len(seqres_df))
    idx += len(seqres_df)

    # ATOM records interleaved with a TER after each chain block
    atom_df = pdb.df["ATOM"].copy()
    ter_rows = []
    for chain_id in atom_df["chain_id"].unique():
        mask = atom_df["chain_id"] == chain_id
        n = int(mask.sum())
        atom_df.loc[mask, "line_idx"] = range(idx, idx + n)
        idx += n
        ter_rows.append({"record_name": "TER",
                         "entry": _ter_entry(atom_df[mask].iloc[-1]),
                         "line_idx": idx})
        idx += 1
    pdb.df["ATOM"] = atom_df

    ter_df = (pd.DataFrame(ter_rows) if ter_rows
              else pd.DataFrame(columns=["record_name", "entry", "line_idx"]))

    hetatm_df = pdb.df["HETATM"].copy()
    hetatm_df["line_idx"] = range(idx, idx + len(hetatm_df))
    idx += len(hetatm_df)
    pdb.df["HETATM"] = hetatm_df

    trailer["line_idx"] = range(idx, idx + len(trailer))

    pdb.df["OTHERS"] = pd.concat([header, seqres_df, ter_df, trailer], ignore_index=True)


def print_info(pdb):
    print("Protein atoms: %d (%d non-hydrogen)" %
          (pdb.df["ATOM"].shape[0], pdb.get('heavy').shape[0]))
    print("Hetero atoms: %d (%d non-hydrogen)" %
          (pdb.df["HETATM"].shape[0],
           pdb.df["HETATM"][pdb.df["HETATM"]["element_symbol"] != "H"].shape[0]))
    print("Coordinate centroid: %f, %f, %f" %
          (pdb.df["ATOM"]["x_coord"].mean(),
           pdb.df["ATOM"]["y_coord"].mean(),
           pdb.df["ATOM"]["z_coord"].mean()))
    print("B-factors (min/max/mean)\n\tProtein: %f / %f / %f\n\tHetero: %f / %f / %f" %
          (pdb.df["ATOM"]["b_factor"].min(), pdb.df["ATOM"]["b_factor"].max(),
           pdb.df["ATOM"]["b_factor"].mean(),
           pdb.df["HETATM"]["b_factor"].min(), pdb.df["HETATM"]["b_factor"].max(),
           pdb.df["HETATM"]["b_factor"].mean()))
    print("Occupancies (min/max/mean)\n\tProtein: %f / %f / %f\n\tHetero: %f / %f / %f" %
          (pdb.df["ATOM"]["occupancy"].min(), pdb.df["ATOM"]["occupancy"].max(),
           pdb.df["ATOM"]["occupancy"].mean(),
           pdb.df["HETATM"]["occupancy"].min(), pdb.df["HETATM"]["occupancy"].max(),
           pdb.df["HETATM"]["occupancy"].mean()))


@click.command()
@click.argument("input", type=click.Path(exists=True))
@click.argument("output", type=click.Path(), required=False, default=None)
@click.option("--renumber-atoms", is_flag=True,
              help="Renumber ATOM and HETATM records from 1 (WILL DELETE CONECT RECORDS)")
@click.option("--sort-residues", is_flag=True,
              help="Sort residues by residue number within each chain, preserving atom order")
@click.option("--renumber-residues", "renumber_residues",
              is_flag=False, flag_value="1", default=None, metavar="N",
              help="Renumber residues per chain from N (default: 1); "
                   "prefix +/- to shift all residue numbers by N instead")
@click.option("--renumber-hetatm-residues", "renumber_hetatm_residues",
              is_flag=False, flag_value="1", default=None, metavar="N",
              help="Renumber HETATM-only chain residues from N (requires --chain); "
                   "same syntax as --renumber-residues")
@click.option("--chain", "chain_opt", metavar="CHAINS",
              help="Comma-separated chain IDs to scope operations to; "
                   "automatically expanded to all chains with identical sequence")
@click.option("--seqres", "seqres_fasta", metavar="FASTA", type=click.Path(exists=True),
              help="Set SEQRES records from a FASTA file (record names must be chain IDs)")
@click.option("--bfactor", metavar="PDB",
              help="Update b-factors from source PDB")
@click.option("--occupancy", metavar="VALUE",
              help="Set all occupancies to a single value")
def main(input, output, renumber_atoms, sort_residues, renumber_residues,
         renumber_hetatm_residues, chain_opt, seqres_fasta, bfactor, occupancy):
    pdb = PandasPdb().read_pdb(input)
    if output is None:
        print_info(pdb)
        return 0
    target_chains = _get_target_chains(pdb, chain_opt)
    if sort_residues:
        pdb.df["ATOM"] = _sort_residues(pdb.df["ATOM"], chain_ids=target_chains)
        pdb.df["HETATM"] = _sort_residues(pdb.df["HETATM"], chain_ids=target_chains)
        pdb.df["ATOM"].loc[:, "atom_number"] = pdb.df["ATOM"].index.values + 1
        pdb.df["HETATM"].loc[:, "atom_number"] = (pdb.df["HETATM"].index.values
                                                  + len(pdb.df["ATOM"]) + 1)
    if renumber_residues is not None:
        if renumber_residues[0] in ("+", "-"):
            shift = int(renumber_residues)
            if target_chains is None:
                pdb.df["ATOM"].loc[:, "residue_number"] += shift
            else:
                mask = pdb.df["ATOM"]["chain_id"].isin(target_chains)
                pdb.df["ATOM"].loc[mask, "residue_number"] += shift
        else:
            start = int(renumber_residues)
            pdb.df["ATOM"] = _renumber_residues(pdb.df["ATOM"], start, chain_ids=target_chains)
        pdb.df["HETATM"] = _resolve_hetatm_collisions(pdb.df["HETATM"], pdb.df["ATOM"])
    if renumber_hetatm_residues is not None:
        if chain_opt is None:
            click.echo("Warning: --renumber-hetatm-residues has no effect without --chain",
                       err=True)
        else:
            atom_chain_ids = set(pdb.df["ATOM"]["chain_id"].unique())
            hetatm_chain_ids = set(pdb.df["HETATM"]["chain_id"].unique())
            target_hetatm_only = (target_chains or set()) & (hetatm_chain_ids - atom_chain_ids)
            if target_hetatm_only:
                if renumber_hetatm_residues[0] in ("+", "-"):
                    shift = int(renumber_hetatm_residues)
                    mask = pdb.df["HETATM"]["chain_id"].isin(target_hetatm_only)
                    pdb.df["HETATM"].loc[mask, "residue_number"] += shift
                else:
                    start = int(renumber_hetatm_residues)
                    pdb.df["HETATM"] = _renumber_residues(pdb.df["HETATM"], start,
                                                          chain_ids=target_hetatm_only)
    if renumber_atoms:
        pdb.df["ATOM"]["atom_number"] = pdb.df["ATOM"].index.values + 1
        pdb.df["HETATM"]["atom_number"] = (pdb.df["HETATM"].index.values
                                           + pdb.df["ATOM"].shape[0])
    if occupancy is not None:
        pdb.df["ATOM"].loc[:, "occupancy"] = float(occupancy)
        if "HETATM" in pdb.df:
            pdb.df["HETATM"].loc[:, "occupancy"] = float(occupancy)
    if bfactor is not None:
        bfact = PandasPdb().read_pdb(bfactor)
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
    if (sort_residues or seqres_fasta or renumber_atoms
            or renumber_residues is not None or renumber_hetatm_residues is not None):
        new_seqres = None
        if seqres_fasta:
            fasta_seqs = list(_parse_fasta(seqres_fasta).values())
            pdb_chains = list(pdb.df["ATOM"]["chain_id"].unique())
            if len(fasta_seqs) != len(pdb_chains):
                click.echo("Warning: %d FASTA sequence(s) but %d chain(s) in PDB"
                           % (len(fasta_seqs), len(pdb_chains)), err=True)
            chains = dict(zip(pdb_chains, fasta_seqs))
            new_seqres = _seqres_entries(chains)
        _rebuild_line_idx(pdb, new_seqres=new_seqres)
    pdb.to_pdb(output)
    return 0


if __name__ == "__main__":
    sys.exit(main(standalone_mode=False))
