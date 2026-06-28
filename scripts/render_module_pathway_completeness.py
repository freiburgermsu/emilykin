#!/usr/bin/env python3
"""Per-phase heatmaps of average pathway completeness across co-occurrence-network
Louvain modules.

Rows    = the 8 pathways whose completeness is summarised in ko_pathway_summary.csv
          (the pathways behind the gene_abundance analysis): Denitrification, DNRA,
          PolyP, Phosphate, PHA, Glycogen, Acetate, Propionate/PHV.
Columns = the Louvain modules defined in the given phase
          (network/network_module_membership_p_value_FDR_phase{N}.json).
Element = mean, over the organisms in that module, of the pathway's completeness
          (genes_present / genes_total).  Each module member iterativeID maps 1:1
          to a MAG via mag_iterativeID_old_to_new.json, and the MAG's per-pathway
          completeness comes from ko_pathway_summary.csv.

Phase II has no co-occurrence modules (its per-phase network had no FDR-passing
edges), so only phases I, III, IV and V are produced.

Output (module_pathway_completeness/):
  module_pathway_completeness_phase{N}.png   heatmap (pathways x modules)
  module_pathway_completeness_phase{N}.csv   the underlying matrix
Run:  ~/Documents/py_venv/bin/python scripts/render_module_pathway_completeness.py
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "module_pathway_completeness"
OUT.mkdir(exist_ok=True)

PHASES = ["I", "III", "IV", "V"]          # phase II: no FDR-passing network -> no modules
PATHWAYS = [                              # (ko_pathway_summary.csv prefix, display label)
    ("Denitrification", "Denitrification"),
    ("Nitrogen_other",  "DNRA"),
    ("PolyP",           "PolyP"),
    ("Phosphate",       "Phosphate"),
    ("PHA",             "PHA"),
    ("Glycogen",        "Glycogen"),
    ("Acetate",         "Acetate"),
    ("Propionate_PHV",  "Propionate/PHV"),
]

# --- per-MAG pathway completeness (genes_present / genes_total) ---
comp = {}
with open(ROOT / "ko_pathway_summary.csv") as f:
    for r in csv.DictReader(f):
        d = {}
        for pre, _ in PATHWAYS:
            tot = float(r[f"{pre}_genes_total"])
            d[pre] = float(r[f"{pre}_genes_present"]) / tot if tot else np.nan
        comp[r["MAG"]] = d

iid2mag = {v: k for k, v in json.load(open(ROOT / "mag_iterativeID_old_to_new.json")).items()}


def module_matrix(modules):
    """modules: list of (key, members[iterativeIDs]); returns (matrix, n_mapped_per_module)."""
    M = np.full((len(PATHWAYS), len(modules)), np.nan)
    nmap = []
    for j, (_k, members) in enumerate(modules):
        mags = [iid2mag[m] for m in members if m in iid2mag]
        nmap.append(len(mags))
        for i, (pre, _) in enumerate(PATHWAYS):
            vals = [comp[m][pre] for m in mags if m in comp and not np.isnan(comp[m][pre])]
            if vals:
                M[i, j] = float(np.mean(vals))
    return M, nmap


for ph in PHASES:
    mm = json.load(open(ROOT / f"network/network_module_membership_p_value_FDR_phase{ph}.json"))
    keys = sorted((k for k in mm if k.startswith("module")), key=lambda k: int(k.split("_")[1]))
    modules = [(k, mm[k]["members"]) for k in keys]
    M, _ = module_matrix(modules)
    labels = [f"M{k.split('_')[1]}\n(n={mm[k]['size']})" for k in keys]

    fig, ax = plt.subplots(figsize=(max(5.0, 0.78 * len(modules) + 2.6), 5.4))
    im = ax.imshow(M, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(modules)), labels=labels, fontsize=9)
    ax.set_yticks(range(len(PATHWAYS)), labels=[d for _, d in PATHWAYS], fontsize=10)
    ax.set_xlabel("Louvain module (n = MAGs)", fontsize=10)
    ax.set_title(f"Average pathway completeness by co-occurrence module — Phase {ph}",
                 fontsize=11, fontweight="bold")
    for i in range(len(PATHWAYS)):
        for j in range(len(modules)):
            v = M[i, j]
            if np.isnan(v):
                ax.text(j, i, "·", ha="center", va="center", color="0.6", fontsize=9)
            else:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if v > 0.55 else "black")
    # white cell separators
    ax.set_xticks(np.arange(-.5, len(modules), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(PATHWAYS), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.02)
    cb.set_label("avg pathway completeness", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f"module_pathway_completeness_phase{ph}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    with open(OUT / f"module_pathway_completeness_phase{ph}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pathway"] + [f"M{k.split('_')[1]}(n={mm[k]['size']})" for k in keys])
        for i, (_, disp) in enumerate(PATHWAYS):
            w.writerow([disp] + ["" if np.isnan(M[i, j]) else f"{M[i, j]:.4f}"
                                 for j in range(len(modules))])
    print(f"phase {ph}: {len(modules)} modules x {len(PATHWAYS)} pathways "
          f"-> module_pathway_completeness_phase{ph}.png/.csv")
