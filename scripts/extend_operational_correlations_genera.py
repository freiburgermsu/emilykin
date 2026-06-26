"""Genus-aggregated counterpart to extend_operational_correlations.py.

Recomputes every (parameter × organism) Spearman correlation against
genus-level abundances (iterativeIDs summed by ``name.split('.')[0]``) and
writes:
    - correlation_dual_qvalue_table_extended_genera.csv
    - correlation_with_time_partial_extended_genera.csv

Methodology mirrors the ASV pipeline:
    * measured-only Spearman + 7-day-gap interpolated Spearman
    * BH-FDR within the entire genus-level result (pooled across all params/genera)
    * confirmed = (q_interp < 0.05) AND (q_meas < 0.05)
    * time-partial Spearman controlling for abundance_day

Original-pipeline parameters are pulled from
``modeling_files/correlations/correlation_with_time_partial.csv``; new env
params come from ``ValueEnviromental2_manual clean.csv``.
"""
import os
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, pearsonr, t as tdist
from statsmodels.stats.multitest import multipletests

ROOT = '/Users/andrewfreiburger/Documents/Research/EmilyKin'
if not os.path.isdir(ROOT):  # fall back to the repo dir this script lives in (e.g. Linux box)
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
# correlation tables live under correlations/ in this repo (modeling_files/correlations/ on the Mac layout)
CORR_DIR = 'correlations' if os.path.isdir('correlations') else 'modeling_files/correlations'

# ---------- genus-aggregated abundances ----------
ab = pd.read_csv('abundances.csv').set_index('sample')
ab.columns = [c.split('.')[0] for c in ab.columns]
ab = ab.T.groupby(level=0).sum().T  # samples × genera (summed across iterativeIDs)
print(f'genus-level abundance matrix: {ab.shape[0]} samples × {ab.shape[1]} genera')

# ---------- date/abundance_day mappings ----------
perf_raw = pd.read_csv('performance_data_with_sample_ids.csv', low_memory=False)
perf_raw['date'] = pd.to_datetime(perf_raw['date'], errors='coerce')
sid_rows = perf_raw[perf_raw['sample_id'].notna() & perf_raw['date'].notna()]
sample_to_day = dict(zip(sid_rows['sample_id'], sid_rows['abundance_day']))
sample_to_date = dict(zip(sid_rows['sample_id'], sid_rows['date']))

# Genus abundances re-keyed by abundance_day (so we can join with perf/env tables)
ab_by_day = ab.copy()
ab_by_day.index = [sample_to_day.get(s, np.nan) for s in ab.index]
ab_by_day = ab_by_day[ab_by_day.index.notna()].groupby(level=0).mean()
ab_by_day.index = ab_by_day.index.astype(float)

# ---------- ORIGINAL-pipeline parameter set ----------
orig = pd.read_csv(f'{CORR_DIR}/correlation_with_time_partial.csv')
ORIG_PARAMS = sorted(orig['parameter'].unique().tolist())
print(f'original-pipeline params to recompute: {len(ORIG_PARAMS)}')

# Extra performance-CSV params needed for the dbRDA variable set but absent from
# the original pipeline (substrate ratios + aerobic time).
EXTRA_PERF_PARAMS = ['COD:N', 'N:P', 'A_time [minutes]']
ORIG_PARAMS = ORIG_PARAMS + [c for c in EXTRA_PERF_PARAMS if c not in ORIG_PARAMS]
# Pull per-sample values for those params from the performance CSV
perf_cols_in_csv = [c for c in ORIG_PARAMS if c in perf_raw.columns]
missing = [c for c in ORIG_PARAMS if c not in perf_raw.columns]
if missing:
    print(f'  [warn] {len(missing)} original params not found as columns in performance CSV: {missing[:5]} ...')
perf_per_day = (perf_raw[perf_raw['sample_id'].notna() & perf_raw['abundance_day'].notna()]
                [['abundance_day'] + perf_cols_in_csv].copy())
perf_per_day = perf_per_day.set_index('abundance_day').apply(pd.to_numeric, errors='coerce').sort_index()
# Per-day measured (collapse duplicates by mean) and 7-day-gap linear interpolation
perf_meas_day = perf_per_day.groupby(level=0).mean()
perf_interp_day = perf_meas_day.interpolate(method='linear', limit=7, limit_area='inside')

# ---------- NEW env-param set from ValueEnviromental2 ----------
env = pd.read_csv('ValueEnviromental2_manual clean.csv')
env['Date_key'] = pd.to_datetime(env['Date_key'], format='%m/%d/%Y', errors='coerce')
env = env[env['Date_key'].notna()].set_index('Date_key')
env = env.drop(columns=[c for c in env.columns if c == 'Date' or c.startswith('Unnamed')], errors='ignore')
env = env.apply(pd.to_numeric, errors='coerce').select_dtypes(include='number').sort_index()
NEW_PARAMS = list(env.columns)
print(f'new env params: {len(NEW_PARAMS)}')

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
env_meas_day = pd.DataFrame(env_meas_rows).T.sort_index()
env_interp_day = pd.DataFrame(env_interp_rows).T.sort_index()


# ---------- correlation helpers ----------
def _spearman_partial(x, y, z):
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


def _correlate_param_set(meas_day_df: pd.DataFrame, interp_day_df: pd.DataFrame, param_list: list[str], label: str):
    """Yield correlation rows for a parameter set against every genus."""
    out = []
    for param in param_list:
        meas_vec = meas_day_df[param].dropna() if param in meas_day_df.columns else pd.Series(dtype=float)
        interp_vec = interp_day_df[param].dropna() if param in interp_day_df.columns else pd.Series(dtype=float)
        for genus in ab_by_day.columns:
            gen_vec = ab_by_day[genus].dropna()
            cm = meas_vec.index.intersection(gen_vec.index)
            rho_m, p_m = (_pair_corr(meas_vec.loc[cm].values, gen_vec.loc[cm].values)
                          if len(cm) >= 3 else (np.nan, np.nan))
            ci = interp_vec.index.intersection(gen_vec.index)
            rho_i, p_i = (_pair_corr(interp_vec.loc[ci].values, gen_vec.loc[ci].values)
                          if len(ci) >= 3 else (np.nan, np.nan))
            if len(ci) >= 4:
                rho_part, p_part = _spearman_partial(
                    interp_vec.loc[ci].values, gen_vec.loc[ci].values,
                    np.asarray(list(ci), dtype=float))
                n_part = len(ci)
            else:
                rho_part, p_part, n_part = np.nan, np.nan, 0
            out.append({
                'parameter': param, 'ASV': genus,
                'rho_interp': rho_i, 'p_interp': p_i, 'n_interp': len(ci),
                'rho_meas': rho_m, 'p_meas': p_m, 'n_meas': len(cm),
                'rho_partial': rho_part, 'p_partial': p_part, 'n_partial': n_part,
            })
    print(f'  {label}: {len(out)} rows computed')
    return out


# ---------- compute ----------
all_rows = []
print('computing original-pipeline params × genera ...')
all_rows += _correlate_param_set(perf_meas_day, perf_interp_day, perf_cols_in_csv, 'original perf')
print('computing new env params × genera ...')
all_rows += _correlate_param_set(env_meas_day, env_interp_day, NEW_PARAMS, 'new env')

df = pd.DataFrame(all_rows)
print(f'total genus-level rows: {len(df)}')

# Pooled BH-FDR over everything (the genus-level set is its own universe)
for pcol, qcol in [('p_meas', 'q_meas'), ('p_interp', 'q_interp'), ('p_partial', 'q_partial')]:
    mask = df[pcol].notna()
    qs = pd.Series(np.nan, index=df.index, dtype=float)
    if mask.any():
        _, q, _, _ = multipletests(df.loc[mask, pcol].values, alpha=0.05, method='fdr_bh')
        qs.loc[mask] = q
    df[qcol] = qs
df['confirmed'] = (df['q_interp'] < 0.05) & (df['q_meas'] < 0.05)
print(f'  confirmed (q_interp<0.05 AND q_meas<0.05): {int(df["confirmed"].sum())}')
print(f'  time-partial passing (q_partial<0.05): {int((df["q_partial"] < 0.05).sum())}')

# Save
dual_cols = ['parameter', 'ASV', 'rho_interp', 'p_interp', 'n_interp',
             'rho_meas', 'p_meas', 'n_meas', 'q_interp', 'q_meas', 'confirmed']
partial_cols = dual_cols + ['rho_partial', 'p_partial', 'n_partial', 'q_partial']
out_dual = f'{CORR_DIR}/correlation_dual_qvalue_table_extended_genera.csv'
out_partial = f'{CORR_DIR}/correlation_with_time_partial_extended_genera.csv'
df[dual_cols].to_csv(out_dual, index=False)
df[partial_cols].to_csv(out_partial, index=False)
print(f'wrote {out_dual}')
print(f'wrote {out_partial}')
