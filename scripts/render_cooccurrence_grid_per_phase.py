"""Per-phase versions of cooccurrence_network_p_value_FDR_*_grid_layout_*_simplified.png.

Reuses render_network from render_cooccurrence_grid.py. For each phase in
phase_day_ranges.json, filters samples (dropping the first 14 days of the phase,
matching the established hygiene convention), recomputes the FDR-passing
co-occurrence graph, swaps in the per-phase mean relative abundances, and
renders with title_slug='p_value_FDR_phase{N}'.
"""
import json
import os
from collections import defaultdict
from itertools import combinations

import networkx as nx
import numpy as np
from pandas import read_csv
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

import render_cooccurrence_grid as rcg

ROOT = '/Users/andrewfreiburger/Documents/Research/EmilyKin'
if not os.path.isdir(ROOT):  # fall back to the repo root (this script lives in scripts/)
    _here = os.path.dirname(os.path.abspath(__file__))
    ROOT = _here if os.path.exists(os.path.join(_here, 'abundances.csv')) else os.path.dirname(_here)
os.chdir(ROOT)

EXCLUDE_DAYS_HEAD = 14  # drop the first N days of each phase
TARGET_MODULES = 15     # aim for ~this many total Louvain communities (modules + singletons)


def find_gamma_for_target(G, target=TARGET_MODULES):
    """Scan Louvain resolution on the positive-rho subgraph of G; return the
    gamma whose total community count (Louvain modules from positive subgraph
    plus isolated singletons added from G) is closest to `target`."""
    pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d['rho'] > 0]
    G_pos = G.edge_subgraph(pos_edges).copy()
    all_nodes = set(G.nodes())
    if G_pos.number_of_edges() == 0:
        return 1.0

    def total_comms(comms):
        covered = set().union(*comms) if comms else set()
        isolates = len(all_nodes - covered)
        return len(comms) + isolates

    coarse = [round(0.1 * i, 2) for i in range(1, 41)]  # 0.1 .. 4.0 step 0.1
    scored = []
    for g in coarse:
        comms = nx.community.louvain_communities(G_pos, weight='weight', resolution=g, seed=42)
        scored.append((g, total_comms(comms)))
    coarse_best = min(scored, key=lambda x: (abs(x[1] - target), x[0]))
    g0 = coarse_best[0]
    fine = [round(g0 + 0.005 * k, 3) for k in range(-9, 10)]
    fine = [g for g in fine if g > 0]
    refined = []
    for g in fine:
        comms = nx.community.louvain_communities(G_pos, weight='weight', resolution=g, seed=42)
        refined.append((g, total_comms(comms)))
    best = min(refined, key=lambda x: (abs(x[1] - target), x[0]))
    return best[0]

sample_days = json.load(open('sample_days.json'))
phase_day_ranges = json.load(open('phase_day_ranges.json'))


def build_phase_graph(abund):
    """Spearman + BH-FDR over co-occurring pairs; return graph of FDR-passing edges."""
    presence = (abund > 0).astype(int)
    cooc = defaultdict(int)
    for s in presence.itertuples(index=False):
        present = [c for c, v in zip(presence.columns, s) if v]
        for p in combinations(sorted(present), 2):
            cooc[p] += 1
    pair_data = []
    for (a, b), cnt in cooc.items():
        rho, pval = spearmanr(abund[a], abund[b])
        if np.isnan(rho):
            continue
        pair_data.append((a, b, rho, pval, cnt))
    if not pair_data:
        return nx.Graph()
    pvals = np.array([t[3] for t in pair_data])
    rejected, _, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    G = nx.Graph()
    for i, (a, b, rho, pval, cnt) in enumerate(pair_data):
        if not rejected[i]:
            continue
        G.add_edge(a, b, weight=abs(rho), rho=rho, pvalue=pval, cooccurrence=cnt)
    return G


abund_full = read_csv('abundances.csv', header=0).set_index('sample')

for phase, span in phase_day_ranges.items():
    lo = span['start'] + EXCLUDE_DAYS_HEAD
    keep = [s for s in abund_full.index if lo <= int(sample_days[s]) <= span['end']]
    ab = abund_full.loc[keep]
    ab = ab.loc[:, (ab.fillna(0) > 0).any(axis=0)]
    print(f"\n=== Phase {phase} (days {lo}-{span['end']}): {ab.shape[0]} samples × {ab.shape[1]} orgs ===")
    if ab.shape[0] < 3 or ab.shape[1] < 2:
        print(f"  Phase {phase}: insufficient data — skipping")
        continue
    G_phase = build_phase_graph(ab)
    print(f"  G: {G_phase.number_of_nodes()} nodes, {G_phase.number_of_edges()} FDR-passing edges")
    if G_phase.number_of_edges() == 0:
        print(f"  Phase {phase}: no FDR-passing edges — skipping")
        continue
    rcg._mean_rel_abund = ab.mean(axis=0)
    gamma_phase = find_gamma_for_target(G_phase, target=TARGET_MODULES)
    print(f"  Phase {phase}: selected gamma={gamma_phase:g} (targeting ~{TARGET_MODULES} total communities)")
    for _simplified in (True, False):
        rcg.render_network(G_phase, f'p_value_FDR_phase{phase}',
                           grid_layout_mode=True, resolution=gamma_phase, simplified=_simplified)
