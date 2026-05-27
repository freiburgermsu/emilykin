"""Time-partial-confirmed-only version of the extended operational correlation
heatmap (the full figure, NOT abundance-filtered).

Reads correlation_with_time_partial_extended.csv and keeps only cells where the
time-partial Spearman q-value is below 0.05 (i.e., the cells that would carry a
lime-green dot in the full version). Saves as
correlation_heatmap_full_with_partial_extended_time.png.
"""
import os, sys
import pandas as pd

GENERA = '--genera' in sys.argv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from _gao_pao_label_helper import color_axis_labels

ROOT = '/Users/andrewfreiburger/Documents/Research/EmilyKin'
os.chdir(ROOT)

# Same env-param list as the extend script
env = pd.read_csv('ValueEnviromental2_manual clean.csv')
env_numeric_cols = env.drop(columns=[c for c in env.columns
                                     if c in ('Date_key', 'Date') or c.startswith('Unnamed')],
                            errors='ignore').apply(pd.to_numeric, errors='coerce').select_dtypes(include='number').columns
NEW_PARAMS = list(env_numeric_cols)

partial_ext = pd.read_csv(
    'modeling_files/correlations/correlation_with_time_partial_extended_genera.csv'
    if GENERA else
    'modeling_files/correlations/correlation_with_time_partial_extended.csv'
)
ab = pd.read_csv('abundances.csv').set_index('sample')

# Time-partial passing rows (these are the lime-dot cells in the full version)
sel = partial_ext[partial_ext['q_partial'].notna() & (partial_ext['q_partial'] < 0.05)].copy()
print(f'time-partial passing rows: {len(sel)}')

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
orig_params = [p for p in rho_mat.index if p not in NEW_PARAMS]
new_in_mat = [p for p in NEW_PARAMS if p in rho_mat.index]
if len(orig_params) > 1:
    orig_params = [orig_params[i] for i in _cluster_order(rho_mat.loc[orig_params], by_axis=0)]
if len(new_in_mat) > 1:
    new_in_mat = [new_in_mat[i] for i in _cluster_order(rho_mat.loc[new_in_mat], by_axis=0)]
rho_mat = rho_mat.reindex(orig_params + new_in_mat)
print(f'final matrix: {rho_mat.shape[0]} params x {rho_mat.shape[1]} ASVs')

fig_w = max(20, 0.16 * rho_mat.shape[1])
fig_h = max(10, 0.45 * rho_mat.shape[0])
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
sns.heatmap(rho_mat, cmap='RdBu_r', norm=norm,
            cbar_kws={'label': 'Spearman ρ', 'shrink': 0.4},
            linewidths=0, ax=ax)

if new_in_mat:
    ax.axhline(len(orig_params), color='black', linewidth=1.5, linestyle='--')

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'Time-partial-confirmed Spearman ρ (extended with ValueEnviromental2): '
             f'{len(rho_mat.index)} parameters × {len(rho_mat.columns)} ASVs '
             f'(q_partial < 0.05 only)', fontsize=12)
ax.tick_params(axis='x', labelsize=5)
ax.tick_params(axis='y', labelsize=8)
plt.setp(ax.get_xticklabels(), rotation=70, ha='right', rotation_mode='anchor')
color_axis_labels(ax, axis='x')
plt.tight_layout()

out = ('modeling_files/correlations/correlation_heatmap_full_with_partial_extended_time'
       f'{"_genera" if GENERA else ""}.png')
plt.savefig(out, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'wrote {out}')
