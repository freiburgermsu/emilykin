#!/usr/bin/env python3
"""Standalone GAO/PAO co-occurrence network figures (aggregate + per-phase).

Spring-layout renderer + drivers extracted verbatim from the (now partly
disabled) data_processing.py notebook cell, so the network/ GAO/PAO figures
pick up the new iterativeID labels.  Reuses render_cooccurrence_grid.py for the
full co-occurrence graph G, the GAO/PAO category set and the abundance frame.

Writes  cooccurrence_network_{GAOs,PAOs}_min0.5pct.png  and
        cooccurrence_network_{GAOs,PAOs}_phase{N}_min0.5pct.png  to CWD (repo root).

Run:  ~/Documents/py_venv/bin/python scripts/render_gao_pao_networks.py
"""
import os, sys, json
from collections import defaultdict
from itertools import combinations

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from pandas import read_csv
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE if os.path.exists(os.path.join(HERE, "abundances.csv")) else os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
import render_cooccurrence_grid as rcg   # builds full G, GAOs_PAOs, _mean_rel_abund

ABT = 0.005                 # 0.5 % mean rel abundance (abund_thresh)
EXCLUDE_DAYS_BELOW = 15     # first 14 days of each phase excluded (notebook constant)

itx = json.load(open("iterativeID_taxonomy.json"))
pcm = json.load(open("Phylum_color_map.json"))
sd  = json.load(open("sample_days.json"))
phr = json.load(open("phase_day_ranges.json"))
GAOs_PAOs = rcg.GAOs_PAOs
LBL = {"GAOs": "green", "Putative GAOs": "mediumseagreen", "PAOs": "blue",
       "Putative PAOs": "cornflowerblue", "Other PHA storing potential+ function": "red"}
inv = {v: k for k, vs in GAOs_PAOs.items() for v in vs}

def gpc(n):
    t = str(n)
    for org, cat in inv.items():
        if org in t:
            return LBL.get(cat)
    return None

g2i = defaultdict(list)
for _id, _t in itx.items():
    g = _t.get("Genus", "")
    if g:
        g2i[g].append(_id)


def render_network_phase(G_sub, mean_rel, suffix, focus=None):
    if G_sub.number_of_edges() == 0:
        print(f"[{suffix}] no edges"); return
    pos = nx.spring_layout(G_sub, seed=42, iterations=100, k=0.3)
    edges = G_sub.edges(data=True)
    rho_values = [d["rho"] for _, _, d in edges]
    edge_widths = [3 * d["weight"] for _, _, d in edges]
    norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = cm.RdBu
    edge_colors = [cmap(norm(r)) for r in rho_values]
    width = 40
    fig, ax = plt.subplots(figsize=(width, 30))
    scale = 5000 * (width / 10) * 2
    cps = np.sqrt
    node_sizes = [scale * cps(mean_rel.get(n, 0)) for n in G_sub.nodes()]
    node_colors = [pcm.get(itx.get(n, {}).get("Phylum", ""), "lightgray") for n in G_sub.nodes()]
    edgecolors = [(gpc(n) or "none") if focus and n in focus else "none" for n in G_sub.nodes()]
    linewidths = [6 if focus and n in focus else 0 for n in G_sub.nodes()]
    nx.draw_networkx_nodes(G_sub, pos, node_size=node_sizes, node_color=node_colors,
                           edgecolors=edgecolors, linewidths=linewidths, alpha=0.9, ax=ax)
    nx.draw_networkx_edges(G_sub, pos, width=edge_widths, edge_color=edge_colors, alpha=0.85, ax=ax)
    max_a = mean_rel.max() if len(mean_rel) else 0
    for n, (x, y) in pos.items():
        ab = mean_rel.get(n, 0)
        fs = 6 + (14 * cps(ab) / cps(max_a) * (width / 10) if max_a > 0 else 0)
        fs = max(5, min(fs, 48))
        ax.text(x, y, str(n), fontsize=fs, color=gpc(n) or "black", fontweight="bold", ha="center", va="center")
    sm = cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("Spearman ρ", fontsize=10 * (width / 10))
    ax.axis("off"); plt.tight_layout()
    out = f"cooccurrence_network_{suffix}.png"
    plt.savefig(out, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fdr_graph(ab):
    presence = (ab.fillna(0) > 0).astype(int)
    cooc = defaultdict(int)
    for s in presence.itertuples(index=False):
        present = [c for c, v in zip(presence.columns, s) if v]
        for p in combinations(sorted(present), 2):
            cooc[p] += 1
    pd_ = []
    for (a, b), cnt in cooc.items():
        r, pv = spearmanr(ab[a], ab[b])
        if np.isnan(r):
            continue
        pd_.append((a, b, r, pv, cnt))
    if not pd_:
        return nx.Graph()
    pvals = np.array([t[3] for t in pd_])
    rej, _, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    G = nx.Graph()
    for i, (a, b, r, pv, cnt) in enumerate(pd_):
        if rej[i]:
            G.add_edge(a, b, weight=abs(r), rho=r, pvalue=pv, cooccurrence=cnt)
    return G


def render_gao_pao(G, mean_rel, phase_suffix):
    for sn in ("GAOs", "PAOs"):
        focus = [i for g in GAOs_PAOs[sn] for i in g2i.get(g, []) if i in G]
        if not focus:
            print(f"[{sn}{phase_suffix}] no focus nodes in G — skipping"); continue
        high = {n for n in G.nodes() if mean_rel.get(n, 0) > ABT}
        keep = set(focus) | high
        es = [(u, v) for u, v in G.edges()
              if (u in set(focus) or v in set(focus)) and u in keep and v in keep]
        Gs = G.edge_subgraph(es).copy()
        render_network_phase(Gs, mean_rel, f"{sn}{phase_suffix}_min{ABT * 100:g}pct", focus=set(focus))


# ---- aggregate (full node-filtered graph from rcg) ----
print("=== aggregate GAOs/PAOs ===")
render_gao_pao(rcg.G, rcg._mean_rel_abund, "")

# ---- per-phase ----
abf = read_csv("abundances.csv", header=0).set_index("sample")
for phase, span in phr.items():
    lo = max(span["start"], EXCLUDE_DAYS_BELOW)
    kp = [s for s in abf.index if lo <= int(sd[s]) <= span["end"]]
    ab = abf.loc[kp]
    ab = ab.loc[:, (ab.fillna(0) > 0).any(axis=0)]
    print(f"=== Phase {phase} (days {lo}-{span['end']}): {ab.shape[0]} samples × {ab.shape[1]} orgs ===")
    if ab.shape[0] < 3 or ab.shape[1] < 2:
        print(f"  phase {phase}: insufficient data — skipping"); continue
    G = fdr_graph(ab)
    if G.number_of_edges() == 0:
        print(f"  phase {phase}: no FDR-passing edges — skipping"); continue
    render_gao_pao(G, ab.mean(axis=0), f"_phase{phase}")
print("DONE")
