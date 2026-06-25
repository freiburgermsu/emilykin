#!/usr/bin/env python3
"""
Part E — reconcile all three clade artifacts to ONE call set (the Part-D tree
placement) and regenerate the figure's claded-RPKM input.

Single source of truth: out/partD_tree_clades.tsv (tree_clade).
Reconciled set: Clade I = {coasm_bin.55, CAN_5_bin.40}; Clade II = the other 11.

Updates (headline clade column set to the tree call; method columns preserved):
  gene_ab_figure/data/nosz_clades.tsv          (+clade_hmm,clade_diamond,tree_clade,tree_support)
  clade_classify/out/nosz_clades_updated.tsv   (clade := tree)
  clade_classify/out/partA_clade_I_II.tsv      (+tree_clade,tree_support; final_clade := tree)
  gene_ab_figure/data/gene_rpkm_per_sample_claded.tsv  (nosZ rows re-split by reconciled clade)
"""
import csv
from pathlib import Path

CC   = Path(__file__).parent
OUT  = CC / "out"
GFD  = CC.parent / "gene_ab_figure" / "data"

# ── source of truth ─────────────────────────────────────────────────────────
tree = {}
with open(OUT / "partD_tree_clades.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        tree[r["mag"]] = r
def t_clade(mag): return tree[mag]["tree_clade"]
def t_supp(mag):  return tree[mag]["placement_support"]
def t_conf(mag):
    # I/II confidence: Clade II membership rests on the 0.84-support clan (high);
    # Clade I rests on the query's sister support.
    if t_clade(mag) == "II": return "high"
    return "high" if (t_supp(mag) and float(t_supp(mag)) >= 0.70) else "low"

recon = {m: t_clade(m) for m in tree}
print("Reconciled call set (tree):  Clade I =",
      [m for m in recon if recon[m] == "I"], " | Clade II count =",
      sum(v == "II" for v in recon.values()))

# ── 1. gene_ab_figure/data/nosz_clades.tsv ──────────────────────────────────
p = GFD / "nosz_clades.tsv"
rows = list(csv.DictReader(open(p), delimiter="\t"))
fld = list(rows[0]) + ["clade_hmm", "clade_diamond", "tree_clade", "tree_support"]
for r in rows:
    m = r["mag"]
    r["clade_hmm"]     = tree[m]["hmm_clade"]
    r["clade_diamond"] = tree[m]["diamond_clade"]
    r["tree_clade"]    = t_clade(m)
    r["tree_support"]  = t_supp(m)
    r["clade"]         = recon[m]            # headline call := tree
    r["confidence"]    = t_conf(m)
with open(p, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fld, delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"  updated {p}")

# ── 2. clade_classify/out/nosz_clades_updated.tsv ───────────────────────────
p = OUT / "nosz_clades_updated.tsv"
rows = list(csv.DictReader(open(p), delimiter="\t"))
fld = list(rows[0])
for r in rows:
    r["clade"] = recon[r["mag"]]
    if "confidence" in r: r["confidence"] = t_conf(r["mag"])
with open(p, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fld, delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"  updated {p}")

# ── 3. clade_classify/out/partA_clade_I_II.tsv ──────────────────────────────
p = OUT / "partA_clade_I_II.tsv"
rows = list(csv.DictReader(open(p), delimiter="\t"))
fld = list(rows[0])
for c in ("tree_clade", "tree_support"):
    if c not in fld: fld.append(c)
for r in rows:
    m = r["mag"]
    r["tree_clade"]   = t_clade(m)
    r["tree_support"] = t_supp(m)
    r["final_clade"]  = recon[m]             # final := tree
with open(p, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fld, delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"  updated {p}")

# ── 4. gene_rpkm_per_sample_claded.tsv — re-split nosZ rows by reconciled clade
p = GFD / "gene_rpkm_per_sample_claded.tsv"
rows = list(csv.DictReader(open(p), delimiter="\t"))
SAMP = [c for c in rows[0] if c.startswith("CAN_")]
# per-MAG nosZ total RPKM = cladeI + cladeII (one is zero) from current file
tot = {}
for r in rows:
    if r["ko_group"] in ("nosZ_cladeI", "nosZ_cladeII"):
        d = tot.setdefault(r["mag"], {s: 0.0 for s in SAMP})
        for s in SAMP:
            d[s] += float(r[s])
flips = []
for r in rows:
    if r["ko_group"] in ("nosZ_cladeI", "nosZ_cladeII"):
        m = r["mag"]; want = "nosZ_clade" + ("I" if recon.get(m) == "I" else "II")
        before = {s: float(r[s]) for s in SAMP}
        for s in SAMP:
            r[s] = f"{(tot[m][s] if r['ko_group'] == want else 0.0):.4f}"
        if any(abs(float(r[s]) - before[s]) > 1e-9 for s in SAMP):
            flips.append((m, r["ko_group"]))
with open(p, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"  updated {p}   (rows changed: {sorted(set(m for m,_ in flips))})")
print("\nDONE — three artifacts + figure RPKM input now carry the single tree call set.")
