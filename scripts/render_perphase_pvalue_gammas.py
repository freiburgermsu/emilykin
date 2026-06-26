#!/usr/bin/env python3
"""Reproduce EVERY committed per-phase p_value_FDR co-occurrence figure with the
new iterativeID labels, at the exact gamma variants present in network/ (the
committed set mixes an auto-gamma sweep, a uniform gamma=1.115 sweep, and a few
manual gammas).  Phase windowing/graph match render_cooccurrence_grid_per_phase.py
(EXCLUDE_DAYS_HEAD=14).  Skips a target if the auto-run already produced it.
Writes to CWD (repo root)."""
import os, sys, json, re, glob
from collections import defaultdict
from itertools import combinations
import numpy as np
import networkx as nx
import matplotlib; matplotlib.use("Agg")
from pandas import read_csv
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE if os.path.exists(os.path.join(HERE, "abundances.csv")) else os.path.dirname(HERE)
os.chdir(ROOT); sys.path.insert(0, HERE)
import render_cooccurrence_grid as rcg

EXCLUDE_DAYS_HEAD = 14
sample_days = json.load(open("sample_days.json"))
phase_day_ranges = json.load(open("phase_day_ranges.json"))

def build_phase_graph(ab):
    presence = (ab.fillna(0) > 0).astype(int)
    cooc = defaultdict(int)
    for s in presence.itertuples(index=False):
        present = [c for c, v in zip(presence.columns, s) if v]
        for p in combinations(sorted(present), 2):
            cooc[p] += 1
    pd_ = []
    for (a, b), cnt in cooc.items():
        rho, pv = spearmanr(ab[a], ab[b])
        if np.isnan(rho):
            continue
        pd_.append((a, b, rho, pv, cnt))
    if not pd_:
        return nx.Graph()
    pvals = np.array([t[3] for t in pd_])
    rej, _, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    G = nx.Graph()
    for i, (a, b, rho, pv, cnt) in enumerate(pd_):
        if rej[i]:
            G.add_edge(a, b, weight=abs(rho), rho=rho, pvalue=pv, cooccurrence=cnt)
    return G

# committed targets: phase -> {(gamma, simplified)}
targets = defaultdict(set)
for f in glob.glob("network/cooccurrence_network_p_value_FDR_phase*_grid_layout_gamma*.png"):
    m = re.search(r"phase([IVX]+)_min0\.5pct_grid_layout_gamma([0-9.]+)(_simplified)?\.png", os.path.basename(f))
    if m:
        targets[m.group(1)].add((float(m.group(2)), bool(m.group(3))))

abund_full = read_csv("abundances.csv", header=0).set_index("sample")
for phase, span in phase_day_ranges.items():
    if phase not in targets:
        continue
    lo = span["start"] + EXCLUDE_DAYS_HEAD
    keep = [s for s in abund_full.index if lo <= int(sample_days[s]) <= span["end"]]
    ab = abund_full.loc[keep]
    ab = ab.loc[:, (ab.fillna(0) > 0).any(axis=0)]
    if ab.shape[0] < 3 or ab.shape[1] < 2:
        print(f"phase {phase}: insufficient"); continue
    G = build_phase_graph(ab)
    if G.number_of_edges() == 0:
        print(f"phase {phase}: no edges"); continue
    rcg._mean_rel_abund = ab.mean(axis=0)
    for gamma, simp in sorted(targets[phase]):
        out = f"cooccurrence_network_p_value_FDR_phase{phase}_min0.5pct_grid_layout_gamma{gamma:g}{'_simplified' if simp else ''}.png"
        if os.path.exists(out):
            print(f"skip (already rendered): {out}"); continue
        print(f"rendering phase {phase} gamma={gamma:g} simplified={simp}")
        rcg.render_network(G, f"p_value_FDR_phase{phase}", grid_layout_mode=True, resolution=gamma, simplified=simp)
print("DONE")
