"""Extend the operational-vs-ASV correlation analysis with the variables in
ValueEnviromental2_manual clean.csv.

- Joins env vars to abundance samples on date via performance_data_with_sample_ids.csv.
- Recomputes Spearman correlations for each (env_var × ASV) pair in three flavors:
    * measured-only (only abundance_days where the env value is directly observed)
    * 7-day gap-capped interpolated (env timeseries interpolated daily, then subset)
    * time-partial (controlling for abundance_day) on the interpolated series
- Applies BH-FDR within the new rows (preserves original q-values verbatim).
- Marks confirmed = (q_meas<0.05 AND q_interp<0.05), same convention as original.
- Appends new rows to correlation_dual_qvalue_table.csv and
  correlation_with_time_partial.csv, written as *_extended.csv.
- Renders correlation_heatmap_full_with_partial_extended.png mirroring the
  original layout (lime-green dot = time-partial q<0.05).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, pearsonr, t as tdist
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from _gao_pao_label_helper import color_axis_labels, order_param_rows

ROOT = '/Users/andrewfreiburger/Documents/Research/EmilyKin'
if not os.path.isdir(ROOT):
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

# ---------------- inputs ----------------
existing_dual = pd.read_csv('correlations/correlation_dual_qvalue_table.csv')
existing_partial = pd.read_csv('correlations/correlation_with_time_partial.csv')

ab = pd.read_csv('abundances.csv').set_index('sample')

perf = pd.read_csv('performance_data_with_sample_ids.csv', low_memory=False)
perf['date'] = pd.to_datetime(perf['date'], errors='coerce')
sid_rows = perf[perf['sample_id'].notna() & perf['date'].notna()]
sample_to_day = dict(zip(sid_rows['sample_id'], sid_rows['abundance_day']))
sample_to_date = dict(zip(sid_rows['sample_id'], sid_rows['date']))

env = pd.read_csv('ValueEnviromental2_manual clean.csv')
env['Date_key'] = pd.to_datetime(env['Date_key'], format='%m/%d/%Y', errors='coerce')
env = env[env['Date_key'].notna()].set_index('Date_key')
env = env.drop(columns=[c for c in env.columns if c == 'Date' or c.startswith('Unnamed')], errors='ignore')
env = env.apply(pd.to_numeric, errors='coerce').select_dtypes(include='number')
env = env.sort_index()
NEW_PARAMS = list(env.columns)
print(f'env params: {len(NEW_PARAMS)} (after coercing to numeric)')

# Per-abundance_day matrix: measured-only and 7-day gap-capped interpolation
env_full_daily = env.resample('D').mean()
env_interp_daily = env_full_daily.interpolate(method='time', limit=7, limit_area='inside')

env_meas_rows, env_interp_rows = {}, {}
for sid, day in sample_to_day.items():
    if sid not in ab.index:
        continue
    date = sample_to_date[sid]
    if date in env_full_daily.index:
        env_meas_rows[day] = env_full_daily.loc[date]
    if date in env_interp_daily.index:
        env_interp_rows[day] = env_interp_daily.loc[date]
env_meas = pd.DataFrame(env_meas_rows).T.sort_index()
env_interp = pd.DataFrame(env_interp_rows).T.sort_index()
env_meas.index.name = env_interp.index.name = 'abundance_day'

# Extra performance-CSV params needed for the dbRDA variable set but absent from
# the original pipeline (substrate ratios + aerobic time). Keyed by abundance_day.
EXTRA_PERF_PARAMS = [c for c in ['COD:N', 'N:P', 'A_time [minutes]'] if c in perf.columns]
_perf_day = (perf[perf['sample_id'].notna() & perf['abundance_day'].notna()]
             [['abundance_day'] + EXTRA_PERF_PARAMS].copy())
_perf_day = _perf_day.set_index('abundance_day').apply(pd.to_numeric, errors='coerce')
perf_meas = _perf_day.groupby(level=0).mean().sort_index()
perf_interp = perf_meas.interpolate(method='linear', limit=7, limit_area='inside')
print(f'extra perf params: {EXTRA_PERF_PARAMS}')

# Abundances re-keyed by abundance_day for joining
ab_by_day = ab.copy()
ab_by_day.index = [sample_to_day.get(s, np.nan) for s in ab.index]
ab_by_day = ab_by_day[ab_by_day.index.notna()]
ab_by_day = ab_by_day.groupby(level=0).mean()  # if duplicates, average


# ---------------- correlation helpers ----------------
def _spearman_partial(x, y, z):
    """Spearman partial correlation of x and y controlling for z; returns (rho, p)."""
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    try:
        rxy = pearsonr(rx, ry)[0]
        rxz = pearsonr(rx, rz)[0]
        ryz = pearsonr(ry, rz)[0]
    except Exception:
        return float('nan'), float('nan')
    denom = np.sqrt(max(0.0, (1 - rxz ** 2) * (1 - ryz ** 2)))
    if denom < 1e-12:
        return float('nan'), float('nan')
    rho_p = (rxy - rxz * ryz) / denom
    n = len(x)
    if n - 3 <= 0:
        return float(rho_p), float('nan')
    tstat = rho_p * np.sqrt((n - 3) / max(1 - rho_p ** 2, 1e-12))
    p = 2 * (1 - tdist.cdf(abs(tstat), df=n - 3))
    return float(rho_p), float(p)


def _pair_corr(x, y):
    try:
        r, p = spearmanr(x, y)
        if np.isnan(r):
            return float('nan'), float('nan')
        return float(r), float(p)
    except Exception:
        return float('nan'), float('nan')


# ---------------- per-pair correlations ----------------
def _rows_for_params(param_list, meas_df, interp_df):
    out = []
    for param in param_list:
        meas_vec = meas_df[param].dropna() if param in meas_df.columns else pd.Series(dtype=float)
        interp_vec = interp_df[param].dropna() if param in interp_df.columns else pd.Series(dtype=float)
        for asv in ab_by_day.columns:
            asv_vec = ab_by_day[asv].dropna()
            cm = meas_vec.index.intersection(asv_vec.index)
            rho_m, p_m = _pair_corr(meas_vec.loc[cm].values, asv_vec.loc[cm].values) if len(cm) >= 3 else (np.nan, np.nan)
            ci = interp_vec.index.intersection(asv_vec.index)
            rho_i, p_i = _pair_corr(interp_vec.loc[ci].values, asv_vec.loc[ci].values) if len(ci) >= 3 else (np.nan, np.nan)
            if len(ci) >= 4:
                rho_p_partial, p_p_partial = _spearman_partial(
                    interp_vec.loc[ci].values, asv_vec.loc[ci].values,
                    np.asarray(list(ci), dtype=float))
                n_p_partial = len(ci)
            else:
                rho_p_partial, p_p_partial, n_p_partial = np.nan, np.nan, 0
            out.append({
                'parameter': param, 'ASV': asv,
                'rho_interp': rho_i, 'p_interp': p_i, 'n_interp': len(ci),
                'rho_meas': rho_m, 'p_meas': p_m, 'n_meas': len(cm),
                'rho_partial': rho_p_partial, 'p_partial': p_p_partial, 'n_partial': n_p_partial,
            })
    return out

rows = _rows_for_params(NEW_PARAMS, env_meas, env_interp)
rows += _rows_for_params(EXTRA_PERF_PARAMS, perf_meas, perf_interp)
new = pd.DataFrame(rows)
print(f'new rows computed: {len(new)} '
      f'({len(NEW_PARAMS)} env + {len(EXTRA_PERF_PARAMS)} extra-perf params × {len(ab_by_day.columns)} ASVs)')

# BH-FDR within the new rows (preserves original q-values verbatim)
for pcol, qcol in [('p_meas', 'q_meas'), ('p_interp', 'q_interp'), ('p_partial', 'q_partial')]:
    mask = new[pcol].notna()
    qs = pd.Series(np.nan, index=new.index, dtype=float)
    if mask.any():
        _, q, _, _ = multipletests(new.loc[mask, pcol].values, alpha=0.05, method='fdr_bh')
        qs.loc[mask] = q
    new[qcol] = qs
new['confirmed'] = (new['q_interp'] < 0.05) & (new['q_meas'] < 0.05)
n_conf_new = int(new['confirmed'].sum())
n_partial_new = int(((new['q_partial'] < 0.05) & new['confirmed']).sum())
print(f'new confirmed: {n_conf_new}, of which time-partial passes: {n_partial_new}')

# ---------------- concat with existing ----------------
dual_cols = ['parameter', 'ASV', 'rho_interp', 'p_interp', 'n_interp',
             'rho_meas', 'p_meas', 'n_meas', 'q_interp', 'q_meas', 'confirmed']
partial_cols = dual_cols + ['rho_partial', 'p_partial', 'n_partial', 'q_partial']
dual_extended = pd.concat([existing_dual, new[dual_cols]], ignore_index=True)
partial_extended = pd.concat([existing_partial, new[partial_cols]], ignore_index=True)
dual_extended.to_csv('correlations/correlation_dual_qvalue_table_extended.csv', index=False)
partial_extended.to_csv('correlations/correlation_with_time_partial_extended.csv', index=False)
print(f'extended tables: dual {dual_extended.shape}, partial {partial_extended.shape}')

# ---------------- render extended heatmap ----------------
conf = partial_extended[partial_extended['confirmed'] == True].copy()
print(f'extended confirmed (orig + new): {len(conf)}')

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

fig_w = max(20, 0.16 * rho_mat.shape[1])
fig_h = max(10, 0.45 * rho_mat.shape[0])
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
sns.heatmap(rho_mat, cmap='RdBu_r', norm=norm,
            cbar_kws={'label': 'Spearman ρ', 'shrink': 0.4}, linewidths=0, ax=ax)

# Lime-green dot overlay for time-partial-confirmed cells
partial_hits = conf[conf['q_partial'].notna() & (conf['q_partial'] < 0.05)]
asv_to_col = {a: i for i, a in enumerate(rho_mat.columns)}
param_to_row = {p: i for i, p in enumerate(rho_mat.index)}
xs, ys = [], []
for _, r in partial_hits.iterrows():
    if r['ASV'] in asv_to_col and r['parameter'] in param_to_row:
        xs.append(asv_to_col[r['ASV']] + 0.5)
        ys.append(param_to_row[r['parameter']] + 0.5)
ax.scatter(xs, ys, s=12, c='lime', marker='o', edgecolors='none', zorder=10)

# Separator line between performance and environmental rows
if 0 < _n_perf < rho_mat.shape[0]:
    ax.axhline(_n_perf, color='black', linewidth=1.5, linestyle='--')

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'Confirmed Spearman ρ (extended with ValueEnviromental2): '
             f'{len(rho_mat.index)} parameters × {len(rho_mat.columns)} ASVs. '
             f'Lime dots = time-partial q<0.05', fontsize=12)
ax.tick_params(axis='x', labelsize=5)
ax.tick_params(axis='y', labelsize=8)
plt.setp(ax.get_xticklabels(), rotation=70, ha='right', rotation_mode='anchor')
color_axis_labels(ax, axis='x')
plt.tight_layout()

out = 'correlations/correlation_heatmap_full_with_partial_extended.png'
plt.savefig(out, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'wrote {out}')
