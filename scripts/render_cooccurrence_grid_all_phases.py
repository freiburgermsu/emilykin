#!/usr/bin/env python3
"""Comprehensive all-phases co-occurrence network + Louvain modules, MAG-filtered.

Built with the SAME methodology as the per-phase networks
(render_cooccurrence_grid_per_phase.py) but over the union of every phase's
cleaned sample set instead of one phase at a time:

  * restrict the node universe to MAG-representative ASV iterativeIDs
    (NETWORK_MAG_ONLY), exactly as the phase-specific networks do;
  * drop the first 14 days of EACH phase (burn-in) and pool all remaining CAN
    samples across phases I-V;
  * recompute the Spearman + BH-FDR (alpha=0.05) co-occurrence graph on that
    pooled matrix;
  * choose the Louvain resolution whose community count (modules + singletons) is
    closest to ~15, then render with the shared grid-layout helper.

This differs from the existing aggregate network (render_cooccurrence_grid.py),
which uses every day (no per-phase burn-in removal), the precomputed
FDR_passing_pairs node set, and a fixed resolution (1.115 -> 3 modules).

Writes (to the repo root; move into network/ alongside the per-phase files):
  network_module_membership_p_value_FDR_allphases.json
  cooccurrence_network_p_value_FDR_allphases_min0.5pct_grid_layout_gamma{g}[_simplified].png

Run from the repo root:
  ~/Documents/py_venv/bin/python scripts/render_cooccurrence_grid_all_phases.py
"""
import os
os.environ.setdefault("NETWORK_MAG_ONLY", "1")   # MAG-only node universe; set before rcg import

import sys
import json
from collections import defaultdict
from itertools import combinations

import networkx as nx
import numpy as np
from pandas import read_csv
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE if os.path.exists(os.path.join(HERE, "abundances.csv")) else os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
import render_cooccurrence_grid as rcg          # provides render_network (+ builds, but does not render, the aggregate G)

EXCLUDE_DAYS_HEAD = 14
TARGET_MODULES = 15
TITLE = "p_value_FDR_allphases"


def build_graph(abund):
    """Spearman + BH-FDR over co-occurring pairs (mirrors the per-phase builder)."""
    presence = (abund > 0).astype(int)
    cooc = defaultdict(int)
    for s in presence.itertuples(index=False):
        present = [c for c, v in zip(presence.columns, s) if v]
        for p in combinations(sorted(present), 2):
            cooc[p] += 1
    pair_data = []
    for (a, b), cnt in cooc.items():
        rho, pval = spearmanr(abund[a], abund[b])
        if not np.isnan(rho):
            pair_data.append((a, b, rho, pval, cnt))
    if not pair_data:
        return nx.Graph()
    rejected, _, _, _ = multipletests([t[3] for t in pair_data], alpha=0.05, method="fdr_bh")
    G = nx.Graph()
    for i, (a, b, rho, pval, cnt) in enumerate(pair_data):
        if rejected[i]:
            G.add_edge(a, b, weight=abs(rho), rho=rho, pvalue=pval, cooccurrence=cnt)
    return G


def find_gamma_for_target(G, target=TARGET_MODULES):
    """Scan Louvain resolution on the positive-rho subgraph; return gamma whose total
    community count (modules + isolated singletons) is closest to target."""
    pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d["rho"] > 0]
    G_pos = G.edge_subgraph(pos_edges).copy()
    all_nodes = set(G.nodes())
    if G_pos.number_of_edges() == 0:
        return 1.0

    def total_comms(comms):
        covered = set().union(*comms) if comms else set()
        return len(comms) + len(all_nodes - covered)

    coarse = [(g, total_comms(nx.community.louvain_communities(
        G_pos, weight="weight", resolution=g, seed=42))) for g in (round(0.1 * i, 2) for i in range(1, 41))]
    g0 = min(coarse, key=lambda x: (abs(x[1] - target), x[0]))[0]
    fine_gs = [round(g0 + 0.005 * k, 3) for k in range(-9, 10) if g0 + 0.005 * k > 0]
    fine = [(g, total_comms(nx.community.louvain_communities(
        G_pos, weight="weight", resolution=g, seed=42))) for g in fine_gs]
    return min(fine, key=lambda x: (abs(x[1] - target), x[0]))[0]


# --- pool every phase's cleaned sample set (drop first 14 days of each phase) ---
sample_days = json.load(open("sample_days.json"))
phase_day_ranges = json.load(open("phase_day_ranges.json"))
abund_full = read_csv("abundances.csv", header=0).set_index("sample")
mag_ids = set(json.load(open("mag_iterativeID_old_to_new.json")).values())
abund_full = abund_full[[c for c in abund_full.columns if c in mag_ids]]
print(f"[MAG-only] {abund_full.shape[1]} MAG-matching ASV columns")

keep = []
for phase, span in phase_day_ranges.items():
    lo = span["start"] + EXCLUDE_DAYS_HEAD
    keep += [s for s in abund_full.index if lo <= int(sample_days[s]) <= span["end"]]
keep = sorted(set(keep))
ab = abund_full.loc[keep]
ab = ab.loc[:, (ab.fillna(0) > 0).any(axis=0)]
print(f"all-phases pooled: {ab.shape[0]} samples x {ab.shape[1]} organisms "
      f"(union of phases I-V, first {EXCLUDE_DAYS_HEAD} days of each dropped)")

G = build_graph(ab)
print(f"G: {G.number_of_nodes()} nodes, {G.number_of_edges()} FDR-passing edges")
if G.number_of_edges() == 0:
    sys.exit("no FDR-passing edges — nothing to render")
gamma = find_gamma_for_target(G, TARGET_MODULES)
print(f"selected gamma={gamma:g} (targeting ~{TARGET_MODULES} communities)")

rcg._mean_rel_abund = ab.mean(axis=0)
for simplified in (True, False):
    rcg.render_network(G, TITLE, grid_layout_mode=True, resolution=gamma, simplified=simplified)
print("DONE")
