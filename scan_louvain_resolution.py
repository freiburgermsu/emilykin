"""Scan Louvain resolution and report module count/sizes.
Reuses the graph G built by render_cooccurrence_grid.py.
"""
import json
from collections import defaultdict
from itertools import combinations

import networkx as nx
import numpy as np
from pandas import DataFrame, read_csv
from scipy.stats import spearmanr

# --- Same graph construction as render_cooccurrence_grid.py ---
_significantly_connected_organisms = [str(x) for x in np.load('FDR_passing_pairs.npy')]
_significantly_connected_organisms.append('Methanobacteriaceae.1')
abundances_3 = read_csv('abundances.csv', header=0).set_index('sample')
abundances_3.drop([c for c in abundances_3.columns if c not in _significantly_connected_organisms], axis=1, inplace=True)

res = spearmanr(abundances_3.values, axis=0, nan_policy='omit')
rho_mat = np.where(np.isnan(np.asarray(res.correlation)), 0.0, np.asarray(res.correlation))
corr_matrix_2 = DataFrame(rho_mat, index=abundances_3.columns, columns=abundances_3.columns)

_presence = (abundances_3 > 0).astype(int)
_cooccurrence = defaultdict(int)
for s in _presence.itertuples(index=False):
    present = [c for c, v in zip(_presence.columns, s) if v]
    for p in combinations(sorted(present), 2):
        _cooccurrence[p] += 1

G = nx.Graph()
for (a, b), cnt in _cooccurrence.items():
    rho = corr_matrix_2.loc[a, b]
    if np.isnan(rho):
        continue
    G.add_edge(a, b, weight=abs(rho), rho=rho, cooccurrence=cnt)

print(f'G: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')
pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d['rho'] > 0]
G_pos = G.edge_subgraph(pos_edges).copy()
print(f'positive-rho subgraph: {G_pos.number_of_nodes()} nodes, {G_pos.number_of_edges()} edges\n')

print(f'{"gamma":>7} | {"#comms":>7} | {"#multi":>7} | {"#singletons":>11} | sizes (top 12)')
print('-' * 80)
for gamma in [1.075, 1.076, 1.077, 1.078, 1.079, 1.080]:
    comms = nx.community.louvain_communities(G_pos, weight='weight', resolution=gamma, seed=42)
    sizes = sorted([len(c) for c in comms], reverse=True)
    n_multi = sum(1 for s in sizes if s >= 2)
    n_single = sum(1 for s in sizes if s == 1)
    print(f'{gamma:>8.4f} | {len(comms):>7} | {n_multi:>7} | {n_single:>11} | {sizes[:12]}')
