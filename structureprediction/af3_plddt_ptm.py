#!/usr/bin/env python

import argparse
import os
import json
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection

def parse_chain_groups(arg, available_chains):
    """
    Parse chain groups.

    Formats:
      - None or 'auto' -> first chain vs all others
      - 'A_B'           -> group1 = ['A'], group2 = ['B']
      - 'A,B_C,D'       -> group1 = ['A','B'], group2 = ['C','D']
      - 'PEPTIDE_CANDIDATERECEPTOR' -> group1 = ['PEPTIDE'], group2 = ['CANDIDATERECEPTOR']
    """
    if arg is None or arg.lower() == "auto":
        if not available_chains:
            raise ValueError("No chains found in structure")
        g1 = [available_chains[0]]
        g2 = available_chains[1:]
        return g1, g2

    if "_" not in arg:
        raise ValueError(f"--chains must contain '_', e.g. 'A_B' or 'A,B_C,D', got '{arg}'")

    g1_str, g2_str = arg.split("_", 1)

    # Always treat comma-separated tokens as full chain IDs
    g1 = [c for c in g1_str.split(",") if c]
    g2 = [c for c in g2_str.split(",") if c]

    if not g1 or not g2:
        raise ValueError(f"--chains '{arg}' did not define non-empty groups")

    for c in g1 + g2:
        if c not in available_chains:
            raise ValueError(
                f"Chain '{c}' from --chains not found in structure; available: {available_chains}"
            )
    return g1, g2

def main():
    parser = argparse.ArgumentParser(
        description="Compute interface pLDDT for AlphaFold3 outputs from CIF + confidences.json."
    )
    parser.add_argument("cif_path", help="Path to AF3 model CIF file (e.g. *_model.cif)")
    parser.add_argument("conf_json", help="Path to AF3 confidences JSON (e.g. *_confidences.json)")
    parser.add_argument(
        "--cutoff",
        type=float,
        default=10.0,
        help="Distance cutoff in Å for interface definition (default: 10.0)",
    )
    parser.add_argument(
        "--chains",
        type=str,
        default="auto",
        help=(
            "Chain groups as 'A_B' (group1 = A, group2 = B), "
            "'A,B_C,D' (group1 = A,B; group2 = C,D), "
            "or 'PEPTIDE_CANDIDATERECEPTOR'. "
            "Default: 'auto' = first chain vs others."
        ),
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default=None,
        help=(
            "Prefix for output file. Final name is <prefix>_af3_plddt_interface.txt. "
            "If not given, uses CIF basename in current directory."
        ),
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default=None,
        help=(
            "Exact output filename (overrides --out-prefix), e.g. "
            "/path/to/my_af3_interface_plddt.txt"
        ),
    )
    args = parser.parse_args()

    # Decide output path
    if args.out_file is not None:
        out_path = args.out_file
    else:
        if args.out_prefix is None:
            base = os.path.basename(args.cif_path)
            if base.lower().endswith(".cif"):
                base = base[:-4]
            out_prefix = base
        else:
            out_prefix = args.out_prefix
        out_path = f"{out_prefix}_af3_plddt_interface.txt"

    # ---- load AF3 confidences (per-atom pLDDT) ----
    with open(args.conf_json) as f:
        conf = json.load(f)

    atom_chain_ids = np.array(conf["atom_chain_ids"])
    atom_plddts = np.array(conf["atom_plddts"])

    # ---- load structure (CIF) ----
    parser_mmcif = MMCIFParser(QUIET=True)
    structure = parser_mmcif.get_structure("model0", args.cif_path)
    model = next(structure.get_models())

    # chains present
    all_chain_ids = [chain.id for chain in model]
    print("Chains in CIF:", all_chain_ids)

    group1_chains, group2_chains = parse_chain_groups(args.chains, all_chain_ids)
    print("Group 1 chains:", group1_chains)
    print("Group 2 chains:", group2_chains)

    # ---- map atoms in CIF to atom_plddts and aggregate to residues ----
    cif_atoms = list(Selection.unfold_entities(model, "A"))
    if len(cif_atoms) != len(atom_plddts):
        raise ValueError(
            f"Number of atoms in CIF ({len(cif_atoms)}) does not match atom_plddts ({len(atom_plddts)})"
        )

    res_list = []  # (chain_id, residue_obj)
    res_plddt_sum = []
    res_plddt_count = []

    last_res = None
    for atom, pld in zip(cif_atoms, atom_plddts):
        res = atom.get_parent()
        chain = res.get_parent()
        key = (chain.id, id(res))
        if last_res is None or key != (last_res[0], id(last_res[1])):
            res_list.append((chain.id, res))
            res_plddt_sum.append(float(pld))
            res_plddt_count.append(1)
            last_res = (chain.id, res)
        else:
            res_plddt_sum[-1] += float(pld)
            res_plddt_count[-1] += 1

    res_plddt = np.array(res_plddt_sum) / np.array(res_plddt_count)
    print("Number of residues (from CIF/atoms):", len(res_list))

    # ---- chain-level residue counts ----
    chain_res_counts = {}
    for chain_id, res in res_list:
        chain_res_counts[chain_id] = chain_res_counts.get(chain_id, 0) + 1
    print("chain_res_counts:", chain_res_counts)

    # ---- Neighbor search for interface residues ----
    all_atoms = cif_atoms
    ns = NeighborSearch(all_atoms)
    chain_atoms = {chain.id: list(chain.get_atoms()) for chain in model}

    res_index_map = {}
    for idx, (c_id, res) in enumerate(res_list):
        res_index_map[(c_id, id(res))] = idx

    interface_res_indices = set()
    cutoff = args.cutoff

    for ci in group1_chains:
        for cj in group2_chains:
            if ci not in chain_atoms or cj not in chain_atoms:
                continue
            for atom in chain_atoms[ci]:
                neighbors = ns.search(atom.coord, cutoff)
                for nb_atom in neighbors:
                    nb_chain = nb_atom.get_parent().get_parent().id
                    if nb_chain != cj:
                        continue
                    res_i = atom.get_parent()
                    res_j = nb_atom.get_parent()
                    idx_i = res_index_map.get((ci, id(res_i)))
                    idx_j = res_index_map.get((cj, id(res_j)))
                    if idx_i is not None:
                        interface_res_indices.add(idx_i)
                    if idx_j is not None:
                        interface_res_indices.add(idx_j)

    interface_res_indices = sorted(interface_res_indices)
    print("Number of interface residues:", len(interface_res_indices))

    # ---- per-chain interface counts and compact labels ----
    iface_counts_per_chain = {}
    compact_by_chain = {}
    three_to_one = {
        "ALA": "A","CYS": "C","ASP": "D","GLU": "E","PHE": "F","GLY": "G","HIS": "H",
        "ILE": "I","LYS": "K","LEU": "L","MET": "M","ASN": "N","PRO": "P","GLN": "Q",
        "ARG": "R","SER": "S","THR": "T","VAL": "V","TRP": "W","TYR": "Y",
    }

    for idx in interface_res_indices:
        chain_id, res = res_list[idx]
        _, resseq, icode = res.get_id()
        resname = res.get_resname()
        one = three_to_one.get(resname.upper(), resname[0])
        label = f"{resseq}{one}"
        iface_counts_per_chain[chain_id] = iface_counts_per_chain.get(chain_id, 0) + 1
        compact_by_chain.setdefault(chain_id, []).append(label)

    print("Interface residues per chain:", iface_counts_per_chain)

    # ---- pLDDT metrics ----
    avg_plddt_all = float(res_plddt.mean())
    if interface_res_indices:
        avg_plddt_interface = float(res_plddt[interface_res_indices].mean())
    else:
        avg_plddt_interface = float("nan")

    # Chain-specific interface pLDDT for balanced metric
    iface_indices_by_chain = {c: [] for c in all_chain_ids}
    for idx in interface_res_indices:
        chain_id, _ = res_list[idx]
        iface_indices_by_chain[chain_id].append(idx)

    avg_plddt_iface_chain = {}
    for c in all_chain_ids:
        idxs = iface_indices_by_chain[c]
        if idxs:
            avg_plddt_iface_chain[c] = float(res_plddt[idxs].mean())
        else:
            avg_plddt_iface_chain[c] = float("nan")

    valid_chain_vals = [v for v in avg_plddt_iface_chain.values() if not np.isnan(v)]
    if valid_chain_vals:
        avg_plddt_interface_balanced = float(np.mean(valid_chain_vals))
    else:
        avg_plddt_interface_balanced = float("nan")

    print(f"Distance cutoff: {cutoff} Å")
    print("Average pLDDT (all residues):", avg_plddt_all)
    print("Average pLDDT (interface residues, size-weighted):", avg_plddt_interface)
    print("Average pLDDT (interface residues, balanced across chains):", avg_plddt_interface_balanced)
    print("Avg interface pLDDT per chain:", avg_plddt_iface_chain)

    # ---- write text output ----
    with open(out_path, "w") as f:
        f.write(f"CIF: {args.cif_path}\n")
        f.write(f"CONFIDENCES_JSON: {args.conf_json}\n")
        f.write(f"Chains_in_model: {all_chain_ids}\n")
        f.write(f"Group1_chains: {group1_chains}\n")
        f.write(f"Group2_chains: {group2_chains}\n")
        f.write(f"chain_res_counts: {chain_res_counts}\n")
        f.write(f"num_residues: {len(res_list)}\n")
        f.write(f"cutoff: {cutoff}\n")
        f.write(f"n_interface_residues: {len(interface_res_indices)}\n")
        f.write(f"interface_residues_per_chain: {iface_counts_per_chain}\n")
        f.write(f"avg_plddt_all: {avg_plddt_all}\n")
        f.write(f"avg_plddt_interface_size_weighted: {avg_plddt_interface}\n")
        f.write(f"avg_plddt_interface_balanced: {avg_plddt_interface_balanced}\n")
        f.write(f"avg_plddt_interface_per_chain: {avg_plddt_iface_chain}\n")
        f.write("\n")

        for chain_id in sorted(compact_by_chain.keys()):
            labels = ", ".join(compact_by_chain[chain_id])
            f.write(
                f"Interacting residues (cutoff {cutoff} Å) chain {chain_id}: {labels}\n"
            )

        f.write("\n# Interface residues (detailed):\n")
        f.write("# chain\tresseq\ticode\tresname\tplddt\n")

        for idx in interface_res_indices:
            chain_id, res = res_list[idx]
            hetflag, resseq, icode = res.get_id()
            resname = res.get_resname()
            res_pld = float(res_plddt[idx])
            f.write(f"{chain_id}\t{resseq}\t{icode}\t{resname}\t{res_pld:.4f}\n")

    print(f"Wrote AF3 interface summary + residues to: {out_path}")

if __name__ == "__main__":
    main()
