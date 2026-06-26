"""1%-max-abundance-filtered version of correlation_heatmap_full_with_partial_extended.png.

Reads correlation_with_time_partial_extended.csv (the output of
extend_operational_correlations.py), filters ASVs to those with max relative
abundance >= 1%, and saves the heatmap as
correlation_heatmap_full_with_partial_extended_1pc.png.
"""
import os, sys
import numpy as np
import pandas as pd

GENERA = '--genera' in sys.argv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from _gao_pao_label_helper import color_axis_labels, order_param_rows

ROOT = '/Users/andrewfreiburger/Documents/Research/EmilyKin'
if not os.path.isdir(ROOT):  # fall back to repo root on the Linux box
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

ABUND_MIN = 0.01

# ---- abundance filter ----
ab = pd.read_csv('abundances.csv').set_index('sample')
if GENERA:
    ab.columns = [c.split('.')[0] for c in ab.columns]
    ab = ab.T.groupby(level=0).sum().T
keep_asvs = set(ab.columns[ab.max(axis=0) >= ABUND_MIN])
print(f'{"genera" if GENERA else "ASVs"} with max rel. abundance >= {ABUND_MIN*100:g}%: {len(keep_asvs)}')

# ---- load extended correlation table (originals + ValueEnviromental2 new params) ----
partial_ext = pd.read_csv(
    'correlations/correlation_with_time_partial_extended_genera.csv'
    if GENERA else
    'correlations/correlation_with_time_partial_extended.csv'
)

# Same ValueEnviromental2 param list as the extend script
env = pd.read_csv('ValueEnviromental2_manual clean.csv')
env_numeric_cols = env.drop(columns=[c for c in env.columns
                                     if c in ('Date_key', 'Date') or c.startswith('Unnamed')],
                            errors='ignore').apply(pd.to_numeric, errors='coerce').select_dtypes(include='number').columns
NEW_PARAMS = list(env_numeric_cols)

conf = partial_ext[(partial_ext['confirmed'] == True) & (partial_ext['ASV'].isin(keep_asvs))].copy()
print(f'confirmed rows after abundance filter: {len(conf)}')

rho_mat = conf.pivot_table(index='parameter', columns='ASV', values='rho_meas', aggfunc='first')
# Hierarchical clustering — columns across all ASVs/genera, rows within each category
from scipy.cluster.hierarchy import linkage, leaves_list


def _cluster_order(mat, by_axis):
    if mat.shape[by_axis] < 2:
        return list(range(mat.shape[by_axis]))
    X = mat.fillna(0).values
    if by_axis == 1:
        X = X.T
    Z = linkage(X, method='average', metric='euclidean')
    return list(leaves_list(Z))


rho_mat = rho_mat.iloc[:, _cluster_order(rho_mat, by_axis=1)]
# Keep only the dbRDA variable set (+ Ax_time); performance rows on top,
# environmental/operational-driver rows below the dashed separator.
rho_mat, _n_perf = order_param_rows(rho_mat)
print(f'final matrix: {rho_mat.shape[0]} params x {rho_mat.shape[1]} ASVs')

# ---- plot ----
fig_w = max(12, 0.30 * rho_mat.shape[1])
fig_h = max(8, 0.45 * rho_mat.shape[0])
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
sns.heatmap(rho_mat, cmap='RdBu_r', norm=norm,
            cbar_kws={'label': 'Spearman ρ', 'shrink': 0.4},
            linewidths=0.4, linecolor='lightgray', ax=ax)

# Lime-green dot overlay for time-partial-confirmed cells
partial_hits = conf[conf['q_partial'].notna() & (conf['q_partial'] < 0.05)]
asv_to_col = {a: i for i, a in enumerate(rho_mat.columns)}
param_to_row = {p: i for i, p in enumerate(rho_mat.index)}
xs, ys = [], []
for _, r in partial_hits.iterrows():
    if r['ASV'] in asv_to_col and r['parameter'] in param_to_row:
        xs.append(asv_to_col[r['ASV']] + 0.5)
        ys.append(param_to_row[r['parameter']] + 0.5)
ax.scatter(xs, ys, s=18, c='lime', marker='o', edgecolors='none', zorder=10)

# Separator between original and new params
if 0 < _n_perf < rho_mat.shape[0]:
    ax.axhline(_n_perf, color='black', linewidth=1.5, linestyle='--')

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'Confirmed Spearman ρ (extended with ValueEnviromental2): '
             f'{len(rho_mat.index)} parameters × {len(rho_mat.columns)} ASVs '
             f'(max rel. abundance ≥ {ABUND_MIN*100:g}%). Lime dots = time-partial q<0.05',
             fontsize=11)
ax.tick_params(axis='x', labelsize=7)
ax.tick_params(axis='y', labelsize=8)
plt.setp(ax.get_xticklabels(), rotation=70, ha='right', rotation_mode='anchor')
color_axis_labels(ax, axis='x')
plt.tight_layout()

out = ('correlations/correlation_heatmap_full_with_partial_extended_1pc'
       f'{"_genera" if GENERA else ""}.png')
plt.savefig(out, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'wrote {out}')
