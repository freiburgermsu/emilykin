"""Filtered version of correlation_heatmap_confirmed_bio.png:
- restricts rows to the 6 biological KPIs the user requested
- filters ASVs to those with max relative abundance >= 1%
- writes correlation_heatmap_confirmed_bio_1pc_filtered.png in the same
  correlations/ directory
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

from _gao_pao_label_helper import color_axis_labels

ROOT = '/Users/andrewfreiburger/Documents/Research/EmilyKin'
if not os.path.isdir(ROOT):  # fall back to repo root on the Linux box
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

ABUND_MIN = 0.01  # 1% max relative abundance threshold

# user-requested KPIs, mapped to the actual parameter strings in the CSV
KPI_DISPLAY = {
    'peakN2O [mg/L]': 'peak N2O [mg/L]',
    'specific denitrification rates [mg NO2–N g−1\xa0VSS−1 h−1]': 'specific denitrification rate',
    'N removal [N%]': 'N removal (%)',
    'N removal (ppm) [N-ppn]': 'N removal (ppm)',
    'P removal (ppm) [P-ppm]': 'P removal (ppm)',
    'P removal [P%]': 'P removal (%)',
}

# --- inputs ---
dual = pd.read_csv(
    'correlations/correlation_dual_qvalue_table_extended_genera.csv'
    if GENERA else
    'correlations/correlation_dual_qvalue_table.csv'
)
ab = pd.read_csv('abundances.csv').set_index('sample')
if GENERA:
    ab.columns = [c.split('.')[0] for c in ab.columns]
    ab = ab.T.groupby(level=0).sum().T
max_ab = ab.max(axis=0)
keep_asvs = set(max_ab[max_ab >= ABUND_MIN].index)
print(f'{"genera" if GENERA else "ASVs"} with max rel. abundance >= {ABUND_MIN*100:g}%: {len(keep_asvs)}')

# --- filter ---
conf = dual[dual['confirmed'] == True].copy()
conf = conf[conf['parameter'].isin(KPI_DISPLAY.keys())]
conf = conf[conf['ASV'].isin(keep_asvs)]
print(f'confirmed rows after KPI + abundance filter: {len(conf)}')
for p in KPI_DISPLAY:
    n_hits = (conf['parameter'] == p).sum()
    print(f'  {KPI_DISPLAY[p]}: {n_hits} ASV hits')

# pivot: rows = KPI, cols = ASV, values = measured-only rho
rho_mat = conf.pivot_table(index='parameter', columns='ASV', values='rho_meas', aggfunc='first')
rho_mat = rho_mat.rename(index=KPI_DISPLAY)
ordered_kpis = [KPI_DISPLAY[p] for p in KPI_DISPLAY if KPI_DISPLAY[p] in rho_mat.index]
rho_mat = rho_mat.reindex(ordered_kpis)

# order ASVs by their max |rho| across the kept KPIs so visually-strong ones cluster
asv_strength = rho_mat.abs().max(axis=0).sort_values(ascending=False)
rho_mat = rho_mat[asv_strength.index]

print(f'final matrix: {rho_mat.shape[0]} KPIs x {rho_mat.shape[1]} ASVs')

# --- plot ---
fig_w = max(8, 0.32 * rho_mat.shape[1])
fig_h = max(3.5, 0.6 * rho_mat.shape[0] + 1)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
sns.heatmap(rho_mat, cmap='RdBu_r', norm=norm,
            cbar_kws={'label': 'Spearman ρ', 'shrink': 0.7},
            linewidths=0.5, linecolor='lightgray',
            ax=ax)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'Confirmed Spearman ρ (interpolated AND measured-only both q<0.05): '
             f'biological KPIs × ASVs (max rel. abundance ≥ {ABUND_MIN*100:g}%)',
             fontsize=10)
ax.tick_params(axis='x', labelsize=8)
ax.tick_params(axis='y', labelsize=10, rotation=0)
plt.setp(ax.get_xticklabels(), rotation=70, ha='right', rotation_mode='anchor')
color_axis_labels(ax, axis='x')
plt.tight_layout()

out = ('correlations/correlation_heatmap_confirmed_bio_1pc_filtered'
       f'{"_genera" if GENERA else ""}.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'wrote {out}')
