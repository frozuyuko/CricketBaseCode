#!/usr/bin/env python3
"""
fix_pairedreceptorinput.py

Repairs AF3 paired-receptor input JSONs under pairedreceptorinput/.
Everything is verified against the sequence; nothing is shifted blindly.

  1. PTM positions   : check each ptmPosition lands on a valid parent residue;
                       auto-correct an unambiguous off-by-one (0-based -> 1-based).
  2. Disulfides      : any existing bondedAtomPair endpoint on a protein "SG" atom
                       must sit on a Cys; same off-by-one self-correction.
  3. GPA2/GPB5 dimer : add single-NAG N-glycan stubs at the given sequons
                       (verifies N-X-S/T, self-corrects +-1). NO inter-chain bond,
                       because the subunits associate non-covalently.
  4. Bursicon dimer  : report Cys positions per chain; add the ONE inter-chain
                       disulfide only if BURS_INTERCHAIN_PAIR is set (see notes).

Dry-run by default. --apply writes in place; the original is copied to <file>.bak once.
Output is tagged OK / FIX / ADD / FLAG so you can grep it.
"""
import argparse, json, os, re, shutil, sys, functools

# stream output immediately instead of block-buffering to the log file
print = functools.partial(print, flush=True)

# --- PTM CCD code -> allowed parent residue(s). Extend with codes you actually use. ---
PTM_PARENTS = {
    "PCA": set("QE"),   # pyroglutamate (usually N-terminal Gln/Glu)
    "TYS": set("Y"),    # sulfotyrosine
    "PTR": set("Y"),    # phosphotyrosine
    "SEP": set("S"),    # phosphoserine
    "TPO": set("T"),    # phosphothreonine
    "HYP": set("P"),    # 4-hydroxyproline
    "HY3": set("P"),
}

# N-glycosylation sequons, as you reported them (1-based; the script confirms / self-corrects).
GLYCO_SITES = {"GPA2": [47, 52], "GPB5": [8, 26]}
# Cys counts from your dump, used only to confirm the subunit-to-chain assignment.
SUBUNIT_CYS = {"GPA2": 10, "GPB5": 5}
SEQUON_ACCEPTORS = set("ST")

# Bursicon inter-chain disulfide. Leave None until you know the pair:
#   - run the dimer unconstrained once, measure inter-chain SG-SG distances,
#     the bonded pair is the one near ~2.05 A (the odd Cys in each chain), OR
#   - take it from an alignment to a bursicon with mapped connectivity.
# Then set, using the CHAIN IDS AS THEY APPEAR IN THE PAIRED FILE (the script
# prints the Cys map to help), e.g.:
#   BURS_INTERCHAIN_PAIR = (("B", 41), ("C", 38))
BURS_INTERCHAIN_PAIR = None


def is_af3_input(d):
    return isinstance(d, dict) and "sequences" in d and d.get("dialect") == "alphafold3"

# AF3 output JSONs we must never parse (huge per-atom arrays / featurized inputs)
OUTPUT_SUFFIXES = ("_data.json", "_confidences.json", "_summary_confidences.json",
                   "_ranking_scores.json", "_ranking_score.json")

def discover(root):
    """Walk the tree WITHOUT following symlinks; keep only candidate input JSONs."""
    found = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        for fn in filenames:
            if not fn.endswith(".json"):
                continue
            if fn.endswith(OUTPUT_SUFFIXES) or fn.startswith(("ranking_scores", "TERMS")):
                continue
            found.append(os.path.join(dirpath, fn))
    return sorted(found)

def looks_like_input(path):
    """Cheap pre-filter: a real AF3 input has '"sequences"' near the top.
    Confidence/score files do not, so we skip them without parsing megabytes."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(8192)
    except OSError:
        return False
    return b'"sequences"' in head

def atomic_write_json(path, data):
    """Back up once, then write via temp + rename so a concurrent reader
    (e.g. a pending af3pred task) never sees a half-written file."""
    bak = path + ".bak"
    if not os.path.exists(bak):
        shutil.copy2(path, bak)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.replace(tmp, path)  # atomic on the same filesystem


def protein_chains(d):
    out = []
    for e in d.get("sequences", []):
        p = e.get("protein")
        if p and p.get("sequence"):
            out.append((p["id"], p["sequence"], p))
    return out

def collect_ids(d):
    ids = set()
    for e in d.get("sequences", []):
        for kind in ("protein", "ligand", "dna", "rna"):
            ent = e.get(kind)
            if not ent:
                continue
            i = ent.get("id")
            if isinstance(i, list):
                ids.update(i)
            elif i is not None:
                ids.add(i)
    return ids

def seq_of(d, chain_id):
    for cid, seq, _ in protein_chains(d):
        if cid == chain_id:
            return seq
    return None

def resolve_position(seq, pos, allowed):
    """Return corrected 1-based pos whose residue is in `allowed`.
    Tries pos (as-1-based) then pos+1 (i.e. value was 0-based). None if neither,
    'AMBIG' if both are valid and differ."""
    cands = []
    for cand in (pos, pos + 1):
        i = cand - 1
        if 0 <= i < len(seq) and seq[i] in allowed:
            cands.append(cand)
    if not cands:
        return None
    if len(cands) == 2 and cands[0] != cands[1]:
        # pos already valid -> keep it (no real shift); else genuinely ambiguous
        return cands[0] if seq[pos - 1] in allowed else "AMBIG"
    return cands[0]

def resolve_sequon(seq, s):
    """Return corrected 1-based Asn position forming N-X-S/T (X!=P). None if none."""
    for cand in (s, s + 1):
        i = cand - 1
        if 0 <= i < len(seq) and seq[i] == "N":
            x, t = i + 1, i + 2
            if t < len(seq) and seq[x] != "P" and seq[t] in SEQUON_ACCEPTORS:
                return cand
    return None

def out_dir_name(name):
    return re.sub(r"[^a-z0-9_-]", "_", (name or "").lower())


# ---------- fixers ----------
def fix_ptms(d, label):
    changed = False
    for cid, seq, p in protein_chains(d):
        for m in (p.get("modifications") or []):
            t, pos = m.get("ptmType"), m.get("ptmPosition")
            if pos is None:
                continue
            allowed = PTM_PARENTS.get(t)
            if allowed is None:
                cur = seq[pos - 1] if 1 <= pos <= len(seq) else "?"
                print(f"  FLAG  {label} chain {cid}: unknown PTM '{t}'@{pos} (on {cur}) - verify manually")
                continue
            r = resolve_position(seq, pos, allowed)
            cur = seq[pos - 1] if 1 <= pos <= len(seq) else "?"
            if r == "AMBIG":
                print(f"  FLAG  {label} chain {cid}: {t}@{pos} ambiguous (both {pos} and {pos+1} are valid)")
            elif r is None:
                print(f"  FLAG  {label} chain {cid}: {t}@{pos} on '{cur}', no valid {sorted(allowed)} at +-1 - verify manually")
            elif r != pos:
                print(f"  FIX   {label} chain {cid}: {t} {pos}('{cur}') -> {r}('{seq[r-1]}')")
                m["ptmPosition"] = r
                changed = True
    return changed

def fix_disulfides(d, label):
    """Self-correct off-by-one on existing protein-SG bond endpoints (disulfides)."""
    changed = False
    for bond in (d.get("bondedAtomPairs") or []):
        for ep in bond:
            if not (isinstance(ep, list) and len(ep) == 3 and ep[2] == "SG"):
                continue
            cid, pos = ep[0], ep[1]
            seq = seq_of(d, cid)
            if seq is None or not isinstance(pos, int):
                continue
            r = resolve_position(seq, pos, set("C"))
            cur = seq[pos - 1] if 1 <= pos <= len(seq) else "?"
            if r in (None, "AMBIG"):
                print(f"  FLAG  {label}: SG bond on chain {cid} res {pos} ('{cur}') not resolvable to Cys - verify")
            elif r != pos:
                print(f"  FIX   {label}: disulfide chain {cid} {pos}('{cur}') -> {r}(Cys)")
                ep[1] = r
                changed = True
    return changed

def assign_subunit(d, subunit):
    """Find the chain that satisfies all of this subunit's sequons (+-1),
    confirmed by Cys count when possible. Returns (chain_id, {site:corrected}) or (None, reason)."""
    sites = GLYCO_SITES[subunit]
    hits = []
    for cid, seq, _ in protein_chains(d):
        corrected = {}
        if all((c := resolve_sequon(seq, s)) is not None and not corrected.update({s: c})
                for s in sites):
            hits.append((cid, seq, corrected, seq.count("C") == SUBUNIT_CYS[subunit]))
    if not hits:
        return None, f"no chain satisfies sequons {sites}"
    if len(hits) == 1:
        return hits[0][0], hits[0][2]
    conf = [h for h in hits if h[3]]
    if len(conf) == 1:
        return conf[0][0], conf[0][2]
    return None, f"ambiguous: chains {[h[0] for h in hits]} all satisfy {sites}"

def add_gpa2gpb5_glycans(d, label):
    if not d.get("bondedAtomPairs"):
        d["bondedAtomPairs"] = []
    used = collect_ids(d)
    n, changed = 1, False
    for subunit in ("GPA2", "GPB5"):
        cid, info = assign_subunit(d, subunit)
        if cid is None:
            print(f"  FLAG  {label}: cannot place {subunit} glycans - {info}")
            continue
        print(f"  NOTE  {label}: {subunit} = chain {cid}")
        for site, asn in info.items():
            gid = f"G{n}"
            while gid in used:
                n += 1; gid = f"G{n}"
            used.add(gid)
            note = "" if asn == site else f" (input {site} -> {asn})"
            d["sequences"].append({"ligand": {"id": gid, "ccdCodes": ["NAG"]}})
            d["bondedAtomPairs"].append([[cid, asn, "ND2"], [gid, 1, "C1"]])
            print(f"  ADD   {label}: NAG {gid} -> chain {cid} Asn{asn}{note}")
            changed = True; n += 1
    return changed

def handle_bursicon(d, label):
    if not d.get("bondedAtomPairs"):
        d["bondedAtomPairs"] = []
    for cid, seq, _ in protein_chains(d):
        cys = [i + 1 for i, a in enumerate(seq) if a == "C"]
        print(f"  NOTE  {label}: chain {cid} {len(cys)} Cys at {cys}")
    if BURS_INTERCHAIN_PAIR is None:
        print(f"  FLAG  {label}: inter-chain disulfide NOT set "
              f"(set BURS_INTERCHAIN_PAIR after an unconstrained run / alignment)")
        return False
    (cA, rA), (cB, rB) = BURS_INTERCHAIN_PAIR
    for c, r in ((cA, rA), (cB, rB)):
        s = seq_of(d, c)
        if s is None or not (1 <= r <= len(s)) or s[r - 1] != "C":
            print(f"  FLAG  {label}: bursicon bond endpoint chain {c} res {r} is not Cys - aborting bond")
            return False
    d["bondedAtomPairs"].append([[cA, rA, "SG"], [cB, rB, "SG"]])
    print(f"  ADD   {label}: inter-chain disulfide {cA}{rA}-{cB}{rB}")
    return True


def kind_of(path, d):
    tag = (os.path.basename(os.path.dirname(path)) + " " + (d.get("name") or "")).lower()
    if "gpa2" in tag and "gpb5" in tag:
        return "gpa2gpb5"
    if "bursicon" in tag:
        return "bursicon"
    return "peptide"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~/finalcricketnrp/af3/pairedreceptorinput"))
    ap.add_argument("--apply", action="store_true", help="write changes in place (else dry-run)")
    args = ap.parse_args()

    files = discover(args.root)
    print(f"Scanning {len(files)} candidate JSON(s) under {args.root}")

    touched = []
    for i, f in enumerate(files, 1):
        if i % 25 == 0:
            print(f"  ... {i}/{len(files)}")
        if not looks_like_input(f):          # skip outputs/confidences without parsing them
            continue
        try:
            d = json.load(open(f))
        except Exception as e:
            print(f"  FLAG  {os.path.relpath(f, args.root)}: cannot parse JSON ({e})")
            continue
        if not is_af3_input(d):
            continue
        label = os.path.relpath(f, args.root)
        k = kind_of(f, d)

        ch = fix_ptms(d, label)
        ch = fix_disulfides(d, label) or ch
        if k == "gpa2gpb5":
            ch = add_gpa2gpb5_glycans(d, label) or ch
        elif k == "bursicon":
            ch = handle_bursicon(d, label) or ch

        if ch:
            touched.append((f, label, out_dir_name(d.get("name"))))
            if args.apply:
                atomic_write_json(f, d)

    print("\n" + "=" * 60)
    if not touched:
        print("No changes needed.")
    else:
        verb = "MODIFIED" if args.apply else "WOULD MODIFY (dry-run)"
        print(f"{verb} {len(touched)} file(s). Clean each output dir and rerun:")
        for _, label, outname in touched:
            print(f"  - {label}    (output dir: {outname})")
    if not args.apply:
        print("\nThis was a dry-run. Re-run with --apply to write (originals -> *.bak).")


if __name__ == "__main__":
    main()
