#!/usr/bin/env python

import argparse
import os
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection


def build_boltz_token_index(cif_path):
    """Replicate Boltz-2 tokenization so CIF residues align with the flat plddt/pae arrays.

    Boltz-2 emits ONE confidence token per polymer residue (standard OR modified, e.g. PCA
    pyroglutamate), and ONE token per atom of a non-polymer ligand (e.g. each NAG glycan atom).
    BioPython collapses a multi-atom ligand into a single residue, so a naive residue<->token 1:1
    mapping breaks whenever a per-atom ligand is present (total_res != len(plddt)). We walk the CIF
    in file (token) order and record, for every polymer residue, the index of its token.

    Polymer vs ligand is decided by label_seq_id: a numeric value marks a polymer residue, '.' marks
    a non-polymer ligand atom (matching how ipsae.py returns None for ligand atoms).

    Returns
        poly_token_idx : {(auth_chain, auth_seq:int, icode): token_index} for polymer residues
        ligand_chains  : set of auth_chain ids that are non-polymer ligands
        n_tokens       : total token count (must equal len(plddt))
        polymer_order  : list of auth_chain ids in first-seen order (polymer chains only)
    """
    field = {}
    fnum = 0
    poly_token_idx = {}
    ligand_chains = set()
    polymer_order = []
    seen_poly_chain = set()
    n_tokens = 0
    cur_key = None
    with open(cif_path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                field[line.strip().split(".", 1)[1]] = fnum
                fnum += 1
                continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            f = line.split()
            label_seq = f[field["label_seq_id"]]
            chain = f[field["auth_asym_id"]] if "auth_asym_id" in field else f[field["label_asym_id"]]
            if label_seq == ".":
                # non-polymer ligand atom -> its own token (Boltz tokenizes ligands per atom)
                ligand_chains.add(chain)
                n_tokens += 1
                cur_key = None
                continue
            auth_seq = f[field["auth_seq_id"]] if "auth_seq_id" in field else label_seq
            icode = f[field["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in field else "?"
            icode = " " if icode in ("?", ".", "") else icode
            key = (chain, int(auth_seq), icode)
            if key != cur_key:
                # first atom of a new polymer residue -> a single token for the whole residue
                poly_token_idx[key] = n_tokens
                n_tokens += 1
                cur_key = key
                if chain not in seen_poly_chain:
                    seen_poly_chain.add(chain)
                    polymer_order.append(chain)
    return poly_token_idx, ligand_chains, n_tokens, polymer_order


def _norm_icode(icode):
    return " " if icode in ("?", ".", "") else icode


def parse_chain_groups(arg, available_chains):
    if arg is None:
        if not available_chains:
            raise ValueError("No chains found in structure")
        g1 = [available_chains[0]]
        g2 = available_chains[1:]
        return g1, g2

    if "_" not in arg:
        raise ValueError(f"--chains must contain '_', e.g. 'A_B' or 'AB_CD', got '{arg}'")

    g1_str, g2_str = arg.split("_", 1)
    g1 = list(g1_str)
    g2 = list(g2_str)

    for c in g1 + g2:
        if c not in available_chains:
            raise ValueError(
                f"Chain '{c}' from --chains not found in structure; available: {available_chains}"
            )
    return g1, g2

def main():
    parser = argparse.ArgumentParser(
        description="Compute interface pLDDT between chain groups from CIF + pLDDT npz (Boltz-2)."
    )
    parser.add_argument("cif_path", help="Path to model CIF file")
    parser.add_argument("plddt_path", help="Path to pLDDT .npz file (with 'plddt' array)")
    parser.add_argument(
        "--cutoff",
        type=float,
        default=10.0,
        help="Distance cutoff in Å for interface definition (default: 10.0)",
    )
    parser.add_argument(
        "--chains",
        type=str,
        default=None,
        help=(
            "Chain groups as 'A_B' (group1 = A, group2 = B) or 'AB_CD' "
            "(group1 = A,B; group2 = C,D). Default: first polymer chain vs other polymer chains. "
            "Non-polymer ligand chains (e.g. NAG glycans) are always excluded."
        ),
    )
    parser.add_argument(
        "--print-chains",
        action="store_true",
        help="Print the polymer chain ids (first-seen order, ligands excluded) and exit. "
             "Used by the scoring pipeline to drive one interface calc per receptor-peptide pair.",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default=None,
        help=(
            "Prefix for output file. Final name is <prefix>_plddt_interface.txt. "
            "If not given, uses CIF basename in current directory."
        ),
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default=None,
        help=(
            "Exact output filename (overrides --out-prefix), e.g. "
            "/path/to/my_interface.txt"
        ),
    )
    args = parser.parse_args()

    # Build the Boltz token map up front (also tells us which chains are ligands).
    poly_token_idx, ligand_chains, n_tokens, polymer_order = build_boltz_token_index(args.cif_path)

    if args.print_chains:
        print(" ".join(polymer_order))
        return

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
        out_path = f"{out_prefix}_plddt_interface.txt"

    # load pLDDT
    data = np.load(args.plddt_path)
    if "plddt" not in data.files:
        raise ValueError(f"'plddt' array not found in {args.plddt_path}; available keys: {data.files}")
    plddt = data["plddt"]

    if n_tokens != len(plddt):
        raise ValueError(
            f"Token count from CIF ({n_tokens}) and pLDDT length ({len(plddt)}) differ!"
        )

    # load structure
    parser_mmcif = MMCIFParser(QUIET=True)
    structure = parser_mmcif.get_structure("model0", args.cif_path)
    model = next(structure.get_models())

    # polymer chains only (exclude per-atom ligand chains such as NAG glycans)
    all_chain_ids = [chain.id for chain in model if chain.id not in ligand_chains]
    print("Polymer chains in model:", all_chain_ids)
    if ligand_chains:
        print("Excluded ligand chains:", sorted(ligand_chains))

    group1_chains, group2_chains = parse_chain_groups(args.chains, all_chain_ids)
    # never let a ligand chain sneak into the groups
    group1_chains = [c for c in group1_chains if c not in ligand_chains]
    group2_chains = [c for c in group2_chains if c not in ligand_chains]
    print("Group 1 chains:", group1_chains)
    print("Group 2 chains:", group2_chains)

    # Build the polymer residue list (token order) and map each residue to its plddt index.
    chain_res_counts = {}
    res_list = []          # (chain_id, residue)
    res_plddt_idx = []     # plddt/token index aligned with res_list
    for chain in model:
        if chain.id in ligand_chains:
            continue
        count = 0
        for res in chain.get_residues():
            _, resseq, icode = res.get_id()
            key = (chain.id, resseq, _norm_icode(icode))
            tok = poly_token_idx.get(key)
            if tok is None:
                # polymer residue with no token (should not happen); skip defensively
                continue
            res_list.append((chain.id, res))
            res_plddt_idx.append(tok)
            count += 1
        chain_res_counts[chain.id] = count

    res_plddt_idx = np.array(res_plddt_idx, dtype=int)

    print("chain_res_counts:", chain_res_counts)
    total_res = sum(chain_res_counts.values())
    print("total polymer residues:", total_res, "plddt tokens:", len(plddt))

    def res_plddt(i):
        return float(plddt[res_plddt_idx[i]])

    # Neighbor search (over polymer atoms only)
    polymer_atoms = []
    chain_atoms = {}
    for chain in model:
        if chain.id in ligand_chains:
            continue
        atoms = list(chain.get_atoms())
        chain_atoms[chain.id] = atoms
        polymer_atoms.extend(atoms)
    ns = NeighborSearch(polymer_atoms)

    # mapping (chain_id, residue_obj) -> flat res_list index
    res_index_map = {}
    for idx, (c_id, res) in enumerate(res_list):
        res_index_map[(c_id, id(res))] = idx

    interface_res_indices = set()
    cutoff = args.cutoff

    # For each chain in group1 vs each chain in group2
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

    # counts per chain at interface and compact labels (e.g. 92F)
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

    # pLDDT metrics
    avg_plddt_all = float(plddt[res_plddt_idx].mean()) if len(res_plddt_idx) else float("nan")
    if interface_res_indices:
        avg_plddt_interface = float(np.mean([res_plddt(i) for i in interface_res_indices]))
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
            avg_plddt_iface_chain[c] = float(np.mean([res_plddt(i) for i in idxs]))
        else:
            avg_plddt_iface_chain[c] = float("nan")

    valid_chain_vals = [v for v in avg_plddt_iface_chain.values() if not np.isnan(v)]
    if valid_chain_vals:
        avg_plddt_interface_balanced = float(np.mean(valid_chain_vals))
    else:
        avg_plddt_interface_balanced = float("nan")

    print(f"Distance cutoff: {cutoff} Å")
    print("Average pLDDT (all polymer residues):", avg_plddt_all)
    print("Average pLDDT (interface residues, size-weighted):", avg_plddt_interface)
    print("Average pLDDT (interface residues, balanced across chains):", avg_plddt_interface_balanced)
    print("Avg interface pLDDT per chain:", avg_plddt_iface_chain)

    # ----- write everything into one txt file -----
    with open(out_path, "w") as f:
        # summary
        f.write(f"CIF: {args.cif_path}\n")
        f.write(f"PLDDT_NPZ: {args.plddt_path}\n")
        f.write(f"Chains_in_model: {all_chain_ids}\n")
        f.write(f"Group1_chains: {group1_chains}\n")
        f.write(f"Group2_chains: {group2_chains}\n")
        f.write(f"Excluded_ligand_chains: {sorted(ligand_chains)}\n")
        f.write(f"chain_res_counts: {chain_res_counts}\n")
        f.write(f"total_residues: {total_res}\n")
        f.write(f"plddt_length: {len(plddt)}\n")
        f.write(f"cutoff: {cutoff}\n")
        f.write(f"n_interface_residues: {len(interface_res_indices)}\n")
        f.write(f"interface_residues_per_chain: {iface_counts_per_chain}\n")
        f.write(f"avg_plddt_all: {avg_plddt_all}\n")
        f.write(f"avg_plddt_interface_size_weighted: {avg_plddt_interface}\n")
        f.write(f"avg_plddt_interface_balanced: {avg_plddt_interface_balanced}\n")
        f.write(f"avg_plddt_interface_per_chain: {avg_plddt_iface_chain}\n")
        f.write("\n")

        # compact per-chain interface lists
        for chain_id in sorted(compact_by_chain.keys()):
            labels = ", ".join(compact_by_chain[chain_id])
            f.write(
                f"Interacting residues (cutoff {cutoff} Å) chain {chain_id}: {labels}\n"
            )

        f.write("\n# Interface residues (detailed):\n")
        f.write("# chain\tresseq\ticode\tresname\tplddt\n")

        # detailed table
        for idx in interface_res_indices:
            chain_id, res = res_list[idx]
            hetflag, resseq, icode = res.get_id()
            resname = res.get_resname()
            f.write(f"{chain_id}\t{resseq}\t{icode}\t{resname}\t{res_plddt(idx):.4f}\n")

    print(f"Wrote interface summary + residues to: {out_path}")

if __name__ == "__main__":
    main()
