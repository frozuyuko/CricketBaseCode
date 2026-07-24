#Please install af3tools 

#!/usr/bin/env python3
"""
Output layout:
    <out_base>/<family>/<peptide>/<receptor_id>__<peptide>.yaml

    python make_boltz_paired_yamls.py \
        --pep_msa_dir  ~/finalcricketnrp/boltz2/msa/peptide \
        --rec_msa_base ~/finalcricketnrp/boltz2/msa/receptor \
        --out_base     ~/finalcricketnrp/boltz2/pairedreceptorinput \
        --ptms         ~/finalcricketnrp/script/peptide_ptms.tsv
"""

import argparse
import csv
import glob
import json
import os
import sys
from pathlib import Path



DIMER_MAP = {
    "gpa2_gpb5_heterodimer": {"B": "gpa2", "C": "gpb5"},
    "bursicon_dimer":        {"B": "bursicon_alpha", "C": "bursicon_beta"}
}


def find_dimer_json(json_dir, stem):
    for f in glob.glob(os.path.join(os.path.expanduser(json_dir), "*_input.json")):
        if os.path.basename(f)[:-len("_input.json")].lower() == stem:
            return f
    return None


def parse_dimer_json(json_path):
    """Return (pca_by_chain, ss_bonds). ss_bonds: [((c1,r1),(c2,r2)), ...]."""
    data = json.load(open(json_path))
    pca, ss = {}, []
    for ent in data.get("sequences", []):
        p = ent.get("protein")
        if not p:
            continue
        cid = p.get("id")
        for m in (p.get("modifications") or []):
            t = (m.get("ptmType") or "")
            if t == "PCA":
                pca.setdefault(cid, []).append(int(m["ptmPosition"]))
            elif t.lower() in ("disulfide_bond", "disulfide", "ssbond"):
                a, b = int(m["ptmPosition"]), int(m["ptmPosition2"])
                ss.append(((cid, min(a, b)), (cid, max(a, b))))
    for pair in (data.get("bondedAtomPairs") or []):
        try:
            (c1, r1, a1), (c2, r2, a2) = pair
        except Exception:
            continue
        if a1 == "SG" and a2 == "SG":
            ss.append(((c1, int(r1)), (c2, int(r2))))
    seen, out = set(), []
    for x, y in ss:
        key = tuple(sorted([x, y]))
        if key not in seen:
            seen.add(key)
            out.append((x, y))
    return pca, out


def build_dimer_yaml(rec_seq, rec_msa, subunits, pca, ss):
    """receptor A + dimer subunits (each its own chain + msa) + SG-SG bonds."""
    L = ["sequences:",
         "  - protein:", "      id: A",
         f"      sequence: {rec_seq}", f"      msa: {rec_msa}"]
    for (cid, seq, msa) in subunits:
        L += ["  - protein:", f"      id: {cid}",
              f"      sequence: {seq}", f"      msa: {msa}"]
        ppos = pca.get(cid, [])
        if ppos:
            L.append("      modifications:")
            for pos in sorted(set(ppos)):
                L += [f"        - position: {pos}", "          ccd: PCA"]
    if ss:
        L.append("constraints:")
        for ((c1, r1), (c2, r2)) in ss:
            L += ["  - bond:",
                  f"      atom1: [{c1}, {r1}, SG]",
                  f"      atom2: [{c2}, {r2}, SG]"]
    return "\n".join(L) + "\n"


def assemble_dimer(stem, pep_dir, json_dir):
    """Return (subunits, pca, ss) for a dimer, or None if a subunit a3m is missing."""
    subunits, seq_by_chain = [], {}
    for cid, sub in DIMER_MAP[stem].items():
        a3m = pep_dir / f"{sub}.a3m"
        if not a3m.exists():
            print(f"  SKIP {stem}: missing subunit a3m {sub}.a3m")
            return None
        _, seq = read_a3m_query(a3m)
        seq_by_chain[cid] = seq
        subunits.append((cid, seq, str(a3m.resolve())))
    pca, ss = {}, []
    jpath = find_dimer_json(json_dir, stem)
    if jpath:
        pca, ss = parse_dimer_json(jpath)
        good = []
        for ((c1, r1), (c2, r2)) in ss:
            s1, s2 = seq_by_chain.get(c1, ""), seq_by_chain.get(c2, "")
            a1 = s1[r1 - 1] if 1 <= r1 <= len(s1) else "?"
            a2 = s2[r2 - 1] if 1 <= r2 <= len(s2) else "?"
            if a1 != "C" or a2 != "C":
                print(f"  WARN {stem}: SS [{c1},{r1}]-[{c2},{r2}] hits {a1}-{a2}, not C-C -> skipped")
            else:
                good.append(((c1, r1), (c2, r2)))
        ss = good
        kind = "inter+intra" if any(a[0] != b[0] for a, b in ss) else "intra"
        print(f"  {stem}: {len(ss)} disulfide(s) [{kind}] from {os.path.basename(jpath)}")
    else:
        print(f"  WARN {stem}: no dimer input JSON found -> no disulfides written")
    return subunits, pca, ss


# --------------------------------------------------------------------------- #
# .a3m parsing                                                                #
# --------------------------------------------------------------------------- #
def read_a3m_query(a3m_path):
    """Return (header, query_sequence) from the FIRST record of an .a3m.
    The query (first record) carries no insertions/gaps, but we uppercase
    and strip '-' / '.' defensively."""
    header, seq = None, []
    with open(a3m_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                    continue
                break  # second record reached
            if header is not None:
                seq.append(line.strip())
    if header is None:
        return None, None
    s = "".join(seq).replace("-", "").replace(".", "").upper()
    return header, s


# --------------------------------------------------------------------------- #
# PTM table                                                                   #
# --------------------------------------------------------------------------- #
def load_ptms(tsv_path):
    """peptide -> dict(sequence, pyroglutamate(int|None),
                       disulfides[(cys_ordinal_i, cys_ordinal_j), ...],
                       amidate(bool))."""
    ptms = {}
    if tsv_path is None:
        return ptms
    cols = ["peptide", "sequence", "pyroglutamate", "disulfides", "amidate_cterm"]
    with open(tsv_path) as fh:
        rows = (ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#"))
        rd = csv.DictReader(rows, fieldnames=cols, delimiter="\t")
        for row in rd:
            name = (row.get("peptide") or "").strip()
            if not name:
                continue
            pg = (row.get("pyroglutamate") or "").strip()
            ds = (row.get("disulfides") or "").strip()
            am = (row.get("amidate_cterm") or "").strip().lower()
            sq = (row.get("sequence") or "").strip().upper()
            pairs = []
            for tok in (t.strip() for t in ds.split(",") if t.strip()):
                a, b = tok.split("-")
                pairs.append((int(a), int(b)))
            ptms[name] = {
                "sequence": sq or None,
                "pyroglutamate": int(pg) if pg else None,
                "disulfides": pairs,  # RESIDUE POSITIONS (1-based), e.g. (7, 14)
                "amidate": am in ("y", "yes", "true", "1"),
            }
    return ptms


# --------------------------------------------------------------------------- #
# helpers                                                                     #
# --------------------------------------------------------------------------- #
def cys_positions(seq):
    """1-based residue positions of every cysteine."""
    return [i + 1 for i, c in enumerate(seq) if c == "C"]


def validate_disulfides(seq, pairs, pep):
    """Return Cys-Cys disulfides as residue positions. Exactly two cysteines ->
    one possible bond, so snap to them and override the table. >2 Cys -> keep the
    table's pairs, dropping any that don't land on C-C."""
    cys = cys_positions(seq)
    if len(cys) == 2:
        snapped = (cys[0], cys[1])
        if pairs and sorted(pairs) != [snapped]:
            print(f"  FIX  {pep}: disulfide {pairs} -> {cys[0]}-{cys[1]} (only two Cys)")
        elif not pairs:
            print(f"  FIX  {pep}: 2 Cys -> disulfide {cys[0]}-{cys[1]} set from sequence")
        return [snapped]
    bonds = []
    for (i, j) in pairs:
        if not (1 <= i <= len(seq) and 1 <= j <= len(seq)):
            print(f"  WARN {pep}: disulfide {i}-{j} out of range (len {len(seq)}) -> skipped")
            continue
        ai, aj = seq[i - 1], seq[j - 1]
        if ai != "C" or aj != "C":
            print(f"  WARN {pep}: disulfide {i}-{j} hits {ai}-{aj}, not C-C -> skipped")
            continue
        bonds.append((i, j))
    if not pairs and len(cys) > 2:
        print(f"  NOTE {pep}: {len(cys)} Cys, no disulfide in table (set pairs manually)")
    return bonds


def build_yaml(rec_seq, rec_msa, pep_seq, pep_msa,
               pyroglu, bonds, amidate):
    """Return the YAML text for one paired complex."""
    L = ["sequences:"]
    # chain A = receptor
    L += ["  - protein:",
          "      id: A",
          f"      sequence: {rec_seq}",
          f"      msa: {rec_msa}"]
    # chain B = peptide
    L += ["  - protein:",
          "      id: B",
          f"      sequence: {pep_seq}",
          f"      msa: {pep_msa}"]
    if pyroglu:
        L += ["      modifications:",
              f"        - position: {pyroglu}",
              "          ccd: PCA"]
    # amidation -> NH2 ligand on chain D (B/C reserved for peptide chains)
    if amidate:
        L += ["  - ligand:",
              "      id: D",
              "      ccd: NH2"]
    # constraints (disulfides first, then the amide bond)
    cons = []
    for (ri, rj) in bonds:
        cons += ["  - bond:",
                 f"      atom1: [B, {ri}, SG]",
                 f"      atom2: [B, {rj}, SG]"]
    if amidate:
        cterm = len(pep_seq)
        cons += ["  - bond:",
                 f"      atom1: [B, {cterm}, C]",
                 "      atom2: [D, 1, N]"]
    if cons:
        L.append("constraints:")
        L += cons
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# main                                                                        #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pep_msa_dir", required=True,
                    help="dir of peptide .a3m (chain B)")
    ap.add_argument("--rec_msa_base", required=True,
                    help="base dir holding <family>/ subdirs of receptor .a3m (chain A)")
    ap.add_argument("--out_base", required=True,
                    help="output base; YAMLs go to <out_base>/<family>/<peptide>/")
    ap.add_argument("--families", nargs="+", default=["rhodopsin", "secretin"])
    ap.add_argument("--ptms", default=None, help="peptide_ptms.tsv")
    ap.add_argument("--pep_json_dir", default="~/finalcricketnrp/json/peptide",
                    help="annotated input JSONs (dimer disulfide connectivity)")
    ap.add_argument("--no_dimers", action="store_true",
                    help="skip the assembled dimer complexes")
    ap.add_argument("--only", nargs="+", default=None,
                    help="restrict to these peptide names (a3m basenames)")
    ap.add_argument("--receptor", default=None,
                    help="restrict to ONE receptor (for array-over-receptors); "
                         "each SLURM task handles one receptor across all peptides")
    ap.add_argument("--skip", nargs="+",
                    default=["gpa2", "gpb5", "bursicon_alpha", "bursicon_beta"],
                    help="peptide names to exclude (default: dimer subunits, "
                         "which are paired as assembled dimers instead)")
    ap.add_argument("--apply", action="store_true",
                    help="actually write files (default is a dry run)")
    a = ap.parse_args()

    pep_dir = Path(a.pep_msa_dir).expanduser()
    rec_base = Path(a.rec_msa_base).expanduser()
    out_base = Path(a.out_base).expanduser()
    ptms = load_ptms(Path(a.ptms).expanduser() if a.ptms else None)

    pep_files = sorted(pep_dir.glob("*.a3m"))
    if a.only:
        keep = set(a.only)
        pep_files = [p for p in pep_files if p.stem in keep]
    if a.skip:
        drop = set(a.skip)
        pep_files = [p for p in pep_files if p.stem not in drop]
    if not pep_files:
        sys.exit(f"No peptide .a3m found in {pep_dir}")

   
    peptides = []  
    print("=== peptide PTM summary ===")
    for pf in pep_files:
        name = pf.stem
        _, seq = read_a3m_query(pf)
        if not seq:
            print(f"WARN {name}: no query sequence parsed -> skipped")
            continue
        ptm = ptms.get(name)
        if ptm is None:
            cys = cys_positions(seq)
            extra = f"; {len(cys)} Cys at {cys} (disulfides NOT set)" if cys else ""
            print(f"  {name:24s} no table entry{extra}")
            peptides.append((name, seq, str(pf.resolve()), None, [], False))
            continue
        if ptm["sequence"] and ptm["sequence"] != seq:
            print(f"  WARN {name}: table sequence != a3m query")
            print(f"        a3m  : {seq}")
            print(f"        table: {ptm['sequence']}")
            print(f"        -> using a3m query; PTM positions follow the a3m")
        pyroglu = ptm["pyroglutamate"]
        if pyroglu:
            aa = seq[pyroglu - 1] if 1 <= pyroglu <= len(seq) else "?"
            if aa not in ("Q", "E", "X"):
                print(f"  WARN {name}: pyroglutamate at pos {pyroglu} but "
                      f"residue is {aa!r} (expected Q or E)")
        if "X" in seq:
            print(f"  WARN {name}: query contains 'X' "
                  f"(masked residue?) -> Boltz wants a clean canonical query")
        bonds = validate_disulfides(seq, ptm["disulfides"], name)
        print(f"  {name:24s} pGlu={pyroglu or '-':>3}  "
              f"disulfides={bonds if bonds else '-'}  "
              f"amide={'yes' if ptm['amidate'] else 'no'}")
        peptides.append((name, seq, str(pf.resolve()),
                         pyroglu, bonds, ptm["amidate"]))

    dimers = []  
    if not a.no_dimers and not a.only:
        print("\n=== dimers ===")
        for stem in DIMER_MAP:
            res = assemble_dimer(stem, pep_dir, a.pep_json_dir)
            if res is not None:
                dimers.append((stem, *res))

    print("\n=== pairing ===")
    total = 0
    for fam in a.families:
        rec_dir = rec_base / fam
        rec_files = sorted(rec_dir.glob("*.a3m"))
        if not rec_files:
            print(f"WARN family {fam}: no receptor .a3m in {rec_dir} -> skipped")
            continue
        recs = []
        for rf in rec_files:
            _, rseq = read_a3m_query(rf)
            if not rseq:
                print(f"WARN receptor {rf.name}: no query sequence -> skipped")
                continue
            recs.append((rf.stem, rseq, str(rf.resolve())))
        if a.receptor:
            recs = [r for r in recs if r[0] == a.receptor]
            if not recs:
                print(f"  NOTE {fam}: receptor {a.receptor!r} not found here -> skipped")
                continue
        for (pep, pseq, pmsa, pyroglu, bonds, amidate) in peptides:
            out_dir = out_base / fam / pep
            if a.apply:
                out_dir.mkdir(parents=True, exist_ok=True)
            for (rec_id, rseq, rmsa) in recs:
                y = build_yaml(rseq, rmsa, pseq, pmsa, pyroglu, bonds, amidate)
                if a.apply:
                    (out_dir / f"{rec_id}__{pep}.yaml").write_text(y)
                total += 1
        for (stem, subunits, pca, ss) in dimers:
            out_dir = out_base / fam / stem
            if a.apply:
                out_dir.mkdir(parents=True, exist_ok=True)
            for (rec_id, rseq, rmsa) in recs:
                y = build_dimer_yaml(rseq, rmsa, subunits, pca, ss)
                if a.apply:
                    (out_dir / f"{rec_id}__{stem}.yaml").write_text(y)
                total += 1
        n_folders = len(peptides) + len(dimers)
        print(f"  {fam:10s} {len(recs):4d} receptors x {n_folders:3d} mature peptides "
              f"= {len(recs) * n_folders:6d} yaml  ({len(peptides)} single + {len(dimers)} dimer folders)")

    verb = "WROTE" if a.apply else "WOULD WRITE (dry run; pass --apply)"
    print(f"\n{verb} {total} yaml files under {out_base}")