#!/usr/bin/env python3
#Require biopython
import io, os, re, json, sys
from Bio import Phylo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(ROOT, "trees", "neuropeptide_gpcr_tree_main.treefile")
REPS_OUT = os.path.join(ROOT, "trees", "main_tree_reps.treefile")
FASTA_IN = os.path.join(ROOT, "input", "combined_165_ref.fasta")
CONF_OUT = os.path.join(ROOT, "results", "neuropeptide_receptors_confident.fasta")
TSV_OUT = os.path.join(ROOT, "results", "neuropeptide_receptors_confident.tsv")
CENT = os.path.join(ROOT, "input", "cent_map.json")

ALRT_MIN, UFB_MIN = 80.0, 95.0


def is_ref(n):
    return n.startswith("REF__")


def famof(n):
    return n.split("__", 1)[1].split("|")[0] if is_ref(n) else None


def parse_sup(s):
    if not s or "/" not in str(s):
        return (None, None)
    a, b = str(s).split("/")[:2]
    try:
        return (float(a), float(b))
    except ValueError:
        return (None, None)


def main():
    raw = open(MAIN).read()
    t = Phylo.read(io.StringIO(raw), "newick")

    reps = None
    if os.path.exists(CENT):
        reps = set(json.load(open(CENT)).get("cent", []))

    def canon(x):
        xl = x.split("__", 1)[1].lower() if "__" in x else x.lower()
        if xl.startswith("gbext"):
            return xl
        m = re.match(r"(?:gbfm|gb0*)(\d+)-r([a-z])", xl)
        return f"g{int(m.group(1))}-r{m.group(2)}" if m else xl

    def is_rep(tip):
        if reps is None:
            return True
        return canon(tip) in reps

    # ---- 1. prune to representatives + references ----
    drop = [tp for tp in t.get_terminals()
            if not is_ref(tp.name) and not is_rep(tp.name)]
    for tp in drop:
        t.prune(tp)
    Phylo.write(t, REPS_OUT, "newick")
    kept = t.get_terminals()
    nref = sum(1 for tp in kept if is_ref(tp.name))
    print(f"pruned tree: {len(kept)} tips ({nref} ref + {len(kept)-nref} cricket reps) "
          f"-> {os.path.relpath(REPS_OUT, ROOT)}")

    # ---- 2. confident set from the FULL tree ----
    t2 = Phylo.read(io.StringIO(raw), "newick")
    conf = []
    for tip in t2.get_terminals():
        if is_ref(tip.name) or not is_rep(tip.name):
            continue
        path = t2.get_path(tip)
        for anc in reversed([t2.root] + path[:-1]):
            refs = [l.name for l in anc.get_terminals() if is_ref(l.name)]
            if refs:
                alrt, ufb = parse_sup(anc.confidence if anc.confidence is not None else anc.name)
                if alrt is not None and ufb is not None and alrt >= ALRT_MIN and ufb >= UFB_MIN:
                    fams = sorted(set(famof(r) for r in refs))
                    conf.append({"tip": tip.name, "families": fams,
                                 "alrt": alrt, "ufboot": ufb})
                break

    seqs, name, buf = {}, None, []
    for ln in open(FASTA_IN):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = ln[1:].split()[0]
            buf = []
        else:
            buf.append(ln)
    if name:
        seqs[name] = "".join(buf)

    seq_by_canon = {canon(h): (h, s) for h, s in seqs.items()}

    import textwrap
    conf.sort(key=lambda r: (r["families"][0], r["tip"]))
    fa, tsv = [], ["receptor_id\ttree_class\tfamily\tSH_aLRT\tUFBoot\tclade_families\tlength_aa"]
    for r in conf:
        ck = canon(r["tip"])
        hit = seq_by_canon.get(ck)
        if not hit:
            continue
        disp = r["tip"].split("__", 1)[1]
        cls = "rhodopsin" if r["tip"].startswith("CRIrho") else "secretin"
        seq = hit[1]
        fa.append(f">{disp} class={cls} family={r['families'][0]} "
                  f"SH-aLRT={r['alrt']:.1f} UFBoot={r['ufboot']:.0f}")
        fa.append("\n".join(textwrap.wrap(seq, 60)))
        tsv.append(f"{disp}\t{cls}\t{r['families'][0]}\t{r['alrt']:.1f}\t"
                   f"{r['ufboot']:.0f}\t{'/'.join(r['families'])}\t{len(seq)}")
    open(CONF_OUT, "w").write("\n".join(fa) + "\n")
    open(TSV_OUT, "w").write("\n".join(tsv) + "\n")