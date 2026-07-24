#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

DIMER_NAMES = {
    "bursicon_dimer":        {"B": "bursicon_alpha", "C": "bursicon_beta"},
    "gpa2_gpb5_heterodimer": {"B": "gpa2",           "C": "gpb5"},
}

def protein_chains(data):
    """Yield (chain_id, sequence, unpairedMsa) for every protein chain."""
    for entity in data.get("sequences", []):
        prot = entity.get("protein")
        if prot is not None:
            yield prot.get("id", ""), prot.get("sequence", ""), prot.get("unpairedMsa", "")

def main(in_dir, out_dir):
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    for jf in sorted(in_dir.glob("*_data.json")):
        base = jf.name[:-len("_data.json")]          # acp_data.json -> acp
        data = json.loads(jf.read_text())
        name_map = DIMER_NAMES.get(base)
        for cid, seq, msa in protein_chains(data):
            if name_map is not None:                 # dimer: use requested per-chain name
                stem = name_map.get(cid)
                if stem is None:
                    print(f"WARN  {jf.name}: unexpected chain id {cid!r}")
                    continue
            else:                                    # single peptide chain
                stem = base

            a3m = msa if msa.strip() else f">{stem}\n{seq}\n"
            out = out_dir / f"{stem}.a3m"
            out.write_text(a3m)
            written += 1
            print(f"{jf.name}  chain {cid} -> {out.name}")
    print(f"\nwrote {written} a3m files to {out_dir}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Extract per-chain unpaired MSAs from AF3 data.json into .a3m")
    ap.add_argument("-i", "--in_dir",  required=True, help="dir of *_data.json files")
    ap.add_argument("-o", "--out_dir", required=True, help="output dir for .a3m files")
    a = ap.parse_args()
    main(a.in_dir, a.out_dir)