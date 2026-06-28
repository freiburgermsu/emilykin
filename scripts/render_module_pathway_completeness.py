#!/usr/bin/env python3
"""Per-phase heatmaps of average pathway completeness across co-occurrence-network
Louvain modules.

Rows    = 4 denitrification steps + 7 other pathways (11 total).  Denitrification is
          shown disaggregated into its four enzymatic steps instead of as a single
          row; the remaining pathways (DNRA, PolyP, Phosphate, PHA, Glycogen, Acetate,
          Propionate/PHV) are taken whole, as summarised in ko_pathway_summary.csv.
Columns = the Louvain modules defined in the given phase
          (network/network_module_membership_p_value_FDR_phase{N}.json).
Element = mean, over the organisms in that module, of the row's completeness
          (genes_present / genes_total).  Each module member iterativeID maps 1:1
          to a MAG via mag_iterativeID_old_to_new.json, and the MAG's per-pathway
          completeness comes from ko_pathway_summary.csv.

Denitrification disaggregation: the single 10-gene "Denitrification" row is split
into its four enzymatic steps (NO3->NO2->NO->N2O->N2).  The four steps partition the
exact same 10 genes (5 + 2 + 2 + 1), so they are a faithful decomposition of the
original row.  The NO3->NO2 (nitrate reductase) step combines the periplasmic (nap:
napA, napB) and membrane-bound (nar: narG/narH/narI) nitrate reductases; note that
narG/narH (K00370/K00371) are shared with the nitrite-oxidoreductase (nxr) of
nitrifiers, so high values there can reflect nitrite oxidisers rather than
denitrifiers.

Phase II has no co-occurrence modules (its per-phase network had no FDR-passing
edges), so only phases I, III, IV and V are produced.

Output (module_pathway_completeness/):
  module_pathway_completeness_phase{N}.png   heatmap (rows x modules)
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

# Denitrification split into its four enzymatic steps.  Each step lists the exact gene
# tokens used in ko_pathway_summary.csv's Denitrification_present_list / _absent_list;
# together they partition the original 10-gene Denitrification row (5 + 2 + 2 + 1).
# Labels read NO3->NO2->NO->N2O->N2 down the rows, i.e. the denitrification cascade.
DENITRIF_STEPS = [                        # (display label, set of gene tokens)
    ("NO₃→NO₂ (nap/nar)",
     {"napA", "napB", "narG/narZ/nxrA", "narH/narY/nxrB", "narI/narV"}),
    ("NO₂→NO (nir)",   {"nirK", "nirS"}),
    ("NO→N₂O (nor)",   {"norB", "norC"}),
    ("N₂O→N₂ (nos)", {"nosZ"}),
]
# Remaining pathways, taken whole: (ko_pathway_summary.csv prefix, display label).
PATHWAYS = [
    ("Nitrogen_other",  "DNRA"),
    ("PolyP",           "PolyP"),
    ("Phosphate",       "Phosphate"),
    ("PHA",             "PHA"),
    ("Glycogen",        "Glycogen"),
    ("Acetate",         "Acetate"),
    ("Propionate_PHV",  "Propionate/PHV"),
]
ROWS = [lab for lab, _ in DENITRIF_STEPS] + [disp for _, disp in PATHWAYS]


def parse_present_genes(cell):
    """'napA(1); napB(1)' -> {'napA', 'napB'};  '' -> set().  Gene tokens may contain
    '/', and each present entry carries a trailing '(count)' that is stripped."""
    genes = set()
    for tok in cell.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        i = tok.rfind("(")
        genes.add(tok[:i].strip() if i != -1 else tok)
    return genes


# --- per-MAG row completeness (genes_present / genes_total) ---
comp = {}
with open(ROOT / "ko_pathway_summary.csv") as f:
    for r in csv.DictReader(f):
        d = {}
        # denitrification steps: from the per-gene present list (total = step size)
        present = parse_present_genes(r["Denitrification_present_list"])
        for label, genes in DENITRIF_STEPS:
            d[label] = len(genes & present) / len(genes)
        # other pathways: from the *_genes_present / *_genes_total columns
        for pre, disp in PATHWAYS:
            tot = float(r[f"{pre}_genes_total"])
            d[disp] = float(r[f"{pre}_genes_present"]) / tot if tot else np.nan
        comp[r["MAG"]] = d

iid2mag = {v: k for k, v in json.load(open(ROOT / "mag_iterativeID_old_to_new.json")).items()}


def module_matrix(modules):
    """modules: list of (key, members[iterativeIDs]); returns (matrix, n_mapped_per_module)."""
    M = np.full((len(ROWS), len(modules)), np.nan)
    nmap = []
    for j, (_k, members) in enumerate(modules):
        mags = [iid2mag[m] for m in members if m in iid2mag]
        nmap.append(len(mags))
        for i, row in enumerate(ROWS):
            vals = [comp[m][row] for m in mags if m in comp and not np.isnan(comp[m][row])]
            if vals:
                M[i, j] = float(np.mean(vals))
    return M, nmap


for ph in PHASES:
    mm = json.load(open(ROOT / f"network/network_module_membership_p_value_FDR_phase{ph}.json"))
    keys = sorted((k for k in mm if k.startswith("module")), key=lambda k: int(k.split("_")[1]))
    modules = [(k, mm[k]["members"]) for k in keys]
    M, _ = module_matrix(modules)
    labels = [f"M{k.split('_')[1]}\n(n={mm[k]['size']})" for k in keys]

    fig, ax = plt.subplots(figsize=(max(5.0, 0.78 * len(modules) + 2.6),
                                    max(5.4, 0.46 * len(ROWS) + 1.8)))
    im = ax.imshow(M, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(modules)), labels=labels, fontsize=9)
    ax.set_yticks(range(len(ROWS)), labels=ROWS, fontsize=9)
    ax.set_xlabel("Louvain module (n = MAGs)", fontsize=10)
    ax.set_title(f"Average pathway completeness by co-occurrence module — Phase {ph}",
                 fontsize=11, fontweight="bold")
    for i in range(len(ROWS)):
        for j in range(len(modules)):
            v = M[i, j]
            if np.isnan(v):
                ax.text(j, i, "·", ha="center", va="center", color="0.6", fontsize=9)
            else:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if v > 0.55 else "black")
    # white cell separators
    ax.set_xticks(np.arange(-.5, len(modules), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(ROWS), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    # heavier rule under the 4 denitrification-step rows, separating them from the
    # whole-pathway rows below
    ax.axhline(len(DENITRIF_STEPS) - 0.5, color="0.15", lw=2.2)
    # bracket label for the denitrification-step group (clear of the long y labels)
    ax.annotate("Denitrification", xy=(0, (len(DENITRIF_STEPS) - 1) / 2),
                xytext=(-11.5, 0), textcoords="offset fontsize",
                xycoords=("axes fraction", "data"), rotation=90,
                ha="center", va="center", fontsize=9, fontweight="bold",
                annotation_clip=False)
    cb = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.02)
    cb.set_label("avg pathway completeness", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f"module_pathway_completeness_phase{ph}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    with open(OUT / f"module_pathway_completeness_phase{ph}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pathway"] + [f"M{k.split('_')[1]}(n={mm[k]['size']})" for k in keys])
        for i, row in enumerate(ROWS):
            w.writerow([row] + ["" if np.isnan(M[i, j]) else f"{M[i, j]:.4f}"
                                for j in range(len(modules))])
    print(f"phase {ph}: {len(modules)} modules x {len(ROWS)} rows "
          f"-> module_pathway_completeness_phase{ph}.png/.csv")
