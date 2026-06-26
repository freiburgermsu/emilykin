"""Clean time-partial-confirmed version of the extended variable set without
any ASV abundance restriction. Same styling as the _1pc_filtered_time figure
but every ASV in the extended table that has at least one time-partial hit
is kept.

Save: correlation_heatmap_confirmed_bio_extended_time.png
"""
import os, sys
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

env = pd.read_csv('ValueEnviromental2_manual clean.csv')
env_numeric_cols = env.drop(columns=[c for c in env.columns
                                     if c in ('Date_key', 'Date') or c.startswith('Unnamed')],
                            errors='ignore').apply(pd.to_numeric, errors='coerce').select_dtypes(include='number').columns
NEW_PARAMS = list(env_numeric_cols)

partial_ext = pd.read_csv(
    'correlations/correlation_with_time_partial_extended_genera.csv'
    if GENERA else
    'correlations/correlation_with_time_partial_extended.csv'
)
sel = partial_ext[(partial_ext['q_partial'].notna()) & (partial_ext['q_partial'] < 0.05)].copy()
print(f'time-partial passing rows (no abundance filter): {len(sel)}')

rho_mat = sel.pivot_table(index='parameter', columns='ASV', values='rho_meas', aggfunc='first')

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
print(f'final matrix: {rho_mat.shape[0]} parameters x {rho_mat.shape[1]} ASVs')

fig_w = max(12, 0.22 * rho_mat.shape[1])
fig_h = max(6, 0.45 * rho_mat.shape[0] + 1)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
sns.heatmap(rho_mat, cmap='RdBu_r', norm=norm,
            cbar_kws={'label': 'Spearman ρ', 'shrink': 0.5},
            linewidths=0.4, linecolor='lightgray', ax=ax)

if 0 < _n_perf < rho_mat.shape[0]:
    ax.axhline(_n_perf, color='black', linewidth=1.5, linestyle='--')

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'Time-partial-confirmed Spearman ρ (extended; q_partial<0.05): '
             f'{len(rho_mat.index)} parameters × {len(rho_mat.columns)} ASVs',
             fontsize=11)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=9, rotation=0)
plt.setp(ax.get_xticklabels(), rotation=70, ha='right', rotation_mode='anchor')
color_axis_labels(ax, axis='x')
plt.tight_layout()

out = ('correlations/correlation_heatmap_confirmed_bio_extended_time'
       f'{"_genera" if GENERA else ""}.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'wrote {out}')
