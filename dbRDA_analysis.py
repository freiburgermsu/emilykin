"""Python translation of dbRDA_codeR.docx (R/vegan ``capscale`` analysis).

Pipeline:
- Long abundance table → filter to phases I–V → top 10 genera per phase (union)
  → sample × genus relative-abundance matrix Y.
- Two predictor sets X come from a workbook with stage-level Performance and
  Environmental sheets, each broadcast to per-sample rows via the sample's stage.
- Distance-based RDA (capscale) of Y vs each predictor set on Bray–Curtis
  dissimilarity, using PCoA followed by linear-constraint extraction (mirrors
  ``vegan::capscale(Y ~ X, distance="bray")``).
- Biplots saved as PNGs: species labels (pink/lavender), constraint arrows
  (blue for performance, dark green for environment), and sample dots colored
  by stage using the same B-id → stage mapping as the R script.

Run:
    /Users/andrewfreiburger/Documents/venv_microbiome_stable/bin/python dbRDA_analysis.py

Inputs are configurable at the top of the file. Defaults look in the EmilyKin
repo root.
"""
from __future__ import annotations
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova as skbio_permanova

# ----- paths ---------------------------------------------------------------
ROOT = Path('/Users/andrewfreiburger/Documents/Research/EmilyKin')
# R: long_file <- file.path(base_dir, "table_rel_full.csv")
LONG_FILE = ROOT / 'table_rel_full.csv'
if not LONG_FILE.exists():  # fall back to the file actually committed to this repo
    LONG_FILE = ROOT / 'table_rel_export.csv'
# R: xlsx_file <- file.path(base_dir, "dbRDAdataset.xlsx") with sheets Performance & Environmental
XLSX_FILE = ROOT / 'dbRDAdataset.xlsx'

OUT_DIR = ROOT / 'graphs' / 'dbRDA'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ----- stage colors & B-id → stage mapping (matches R) ---------------------
STAGE_COLORS = {
    'I':   '#EE9A49',
    'II':  '#BFD02F',
    'III': '#90BE66',
    'IV':  '#425F28',
    'V':   '#8482CA',
}

# Same GAO/PAO functional categories used by the heatmaps' axis-label coloring
GAOs_PAOs = {
    'PAOs': ['Ca_Accumulibacter', 'Tetrasphaera', 'Dechloromonas', 'Microlunatus',
             'Azonexus', 'Ca_Phosphoribacter'],
    'GAOs': ['Ca_Competibacter', 'Defluviicoccus', 'Propionivibrio', 'Ca_Contendobacter'],
    'Putative PAOs': ['Ca_Obscuribacter', 'Thauera', 'Zoogloea', 'Paracoccus'],
    'Putative GAOs': ['Micropruina', 'Amaricoccus', 'Ca_Glycocaulis', 'Thauera'],
    'Other PHA storing potential+ function': ['Pseudomonas', 'Bacillus', 'Acinetobacter',
                                              'Rhodocyclaceae'],
}
GAO_PAO_LABEL_COLORS = {
    'GAOs': 'green',
    'Putative GAOs': 'mediumseagreen',
    'PAOs': 'blue',
    'Putative PAOs': 'cornflowerblue',
    'Other PHA storing potential+ function': 'red',
}
SPECIES_DEFAULT_COLOR = 'black'
_INVERTED_GAOs_PAOs = {v: k for k, vs in GAOs_PAOs.items() for v in vs}


def species_label_color(genus: str) -> str:
    """Mirror the heatmaps' substring-match label coloring."""
    for org, cat in _INVERTED_GAOs_PAOs.items():
        if org in str(genus):
            return GAO_PAO_LABEL_COLORS.get(cat, SPECIES_DEFAULT_COLOR)
    return SPECIES_DEFAULT_COLOR

def species_category(genus: str) -> str:
    """GAO/PAO functional category for a genus (substring match), or '' if none."""
    for org, cat in _INVERTED_GAOs_PAOs.items():
        if org in str(genus):
            return cat
    return ''

def stage_from_id(sample_id: str) -> str | None:
    """B17–B25 → I, B26–B30 → II, B31–B41 → III, B42–B52 → IV, B53–B60 → V."""
    s = str(sample_id)
    if not s.startswith('B'):
        return None
    try:
        n = int(s[1:])
    except ValueError:
        return None
    if 17 <= n <= 25: return 'I'
    if 26 <= n <= 30: return 'II'
    if 31 <= n <= 41: return 'III'
    if 42 <= n <= 52: return 'IV'
    if 53 <= n <= 60: return 'V'
    return None


# ----- abundance matrix Y --------------------------------------------------
def build_abundance_matrix(long_file: Path) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Return (Y, sample_ids, stages). Mirrors the first R chunk."""
    tab = pd.read_csv(long_file, low_memory=False)
    # filter to phases I–V
    tab_stage = tab[tab['Phase'].isin(['I', 'II', 'III', 'IV', 'V'])].copy()

    # top 10 genera per stage (union)
    top_genus = (
        tab_stage.dropna(subset=['Genus'])
        .groupby(['Phase', 'Genus'], as_index=False)
        .agg(mean_rel_ab=('rel_ab', 'mean'))
    )
    top_per_stage = (
        top_genus.sort_values(['Phase', 'mean_rel_ab'], ascending=[True, False])
        .groupby('Phase', as_index=False)
        .head(10)
    )
    keep_genera = sorted(top_per_stage['Genus'].unique())

    # sample × genus matrix
    abund_long = (
        tab_stage[tab_stage['Genus'].isin(keep_genera)]
        .groupby(['sample', 'Phase', 'Genus'], as_index=False)
        .agg(rel_ab=('rel_ab', 'sum'))
    )
    Y_wide = abund_long.pivot_table(
        index=['sample', 'Phase'], columns='Genus', values='rel_ab', fill_value=0,
    ).reset_index()

    # sort: by phase order (I…V), then sample
    phase_order = {p: i for i, p in enumerate(['I', 'II', 'III', 'IV', 'V'])}
    Y_wide['_phase_ord'] = Y_wide['Phase'].map(phase_order)
    Y_wide = Y_wide.sort_values(['_phase_ord', 'sample']).drop(columns=['_phase_ord'])

    stages = Y_wide['Phase'].tolist()
    sample_ids = Y_wide['sample'].tolist()
    Y = Y_wide.drop(columns=['sample', 'Phase']).reset_index(drop=True)
    Y.index = sample_ids
    return Y, sample_ids, stages


# ----- db-RDA (capscale equivalent) ----------------------------------------
def dbrda(Y: pd.DataFrame, X: pd.DataFrame, metric: str = 'braycurtis') -> dict:
    """Distance-based RDA ≅ vegan::capscale(Y ~ X, distance=<metric>).

    Steps:
      1. Bray–Curtis distance on Y rows.
      2. PCoA → principal coordinates F (positive eigenvalues only).
      3. Center X; OLS-regress F on X → fitted F_hat.
      4. SVD on F_hat → constrained eigenvalues (S²) and axes.
      5. Site / species / biplot scores in vegan ``scaling = 2`` convention.
    """
    assert Y.index.equals(X.index), 'Y and X must share index'
    n = Y.shape[0]

    # 1. Bray–Curtis on samples
    D = beta_diversity(metric, Y.values, ids=list(Y.index))

    # 2. PCoA — keep positive eigenvalues
    pres = pcoa(D, number_of_dimensions=min(n - 1, Y.shape[1] + 10))
    eig = pres.eigvals.values
    F_all = pres.samples.values  # rows = samples, cols = axes; already F = U * sqrt(eig)
    pos = eig > 1e-12
    F = F_all[:, pos]
    eig_pos = eig[pos]
    total_inertia = eig_pos.sum()

    # 3. Center X
    X_vals = X.values.astype(float)
    X_c = X_vals - X_vals.mean(axis=0, keepdims=True)

    # OLS: F = X_c B + E  →  B = (X_c'X_c)^-1 X_c'F
    B, *_ = np.linalg.lstsq(X_c, F, rcond=None)
    F_hat = X_c @ B   # fitted = constrained component

    # 4. SVD on constrained component
    U, S, Vt = np.linalg.svd(F_hat, full_matrices=False)
    k = int((S > 1e-12).sum())
    U, S, Vt = U[:, :k], S[:k], Vt[:k, :]
    eig_c = S ** 2
    constrained_inertia = eig_c.sum()

    # 5. Scores
    # Linear-combination site scores (constrained space): F_hat @ V
    V = Vt.T
    site_lc = F_hat @ V                            # n × k
    site_wa = F @ V                                # n × k (weighted-average)

    # Species scores: correlation of each species column with linear-combination site scores
    Y_c = Y.values - Y.values.mean(axis=0, keepdims=True)
    sd_sp = Y_c.std(axis=0, ddof=0)
    sd_site = site_lc.std(axis=0, ddof=0)
    sd_sp[sd_sp == 0] = 1.0
    sd_site[sd_site == 0] = 1.0
    sp_scores = (Y_c.T @ site_lc) / (n * np.outer(sd_sp, sd_site))

    # Biplot scores (correlation of each X variable with constrained site scores)
    Xc_std = X_c.std(axis=0, ddof=0)
    Xc_std[Xc_std == 0] = 1.0
    bp_scores = (X_c.T @ site_lc) / (n * np.outer(Xc_std, sd_site))

    # vegan scaling = 2 ("species" scaling): species scaled by sqrt(eig / total), sites unscaled
    scale_eig = np.sqrt(eig_c / total_inertia) if total_inertia > 0 else np.ones_like(eig_c)
    sp_scaled = sp_scores * scale_eig
    bp_scaled = bp_scores * scale_eig

    return {
        'sites': pd.DataFrame(site_lc, index=Y.index, columns=[f'CAP{i+1}' for i in range(k)]),
        'sites_wa': pd.DataFrame(site_wa, index=Y.index, columns=[f'CAP{i+1}' for i in range(k)]),
        'species': pd.DataFrame(sp_scaled, index=Y.columns, columns=[f'CAP{i+1}' for i in range(k)]),
        'biplot': pd.DataFrame(bp_scaled, index=X.columns, columns=[f'CAP{i+1}' for i in range(k)]),
        'eig_constrained': pd.Series(eig_c, index=[f'CAP{i+1}' for i in range(k)]),
        'total_inertia': total_inertia,
        'constrained_inertia': constrained_inertia,
        'prop_constrained': constrained_inertia / total_inertia if total_inertia > 0 else np.nan,
    }


# ----- PERMANOVA companions ------------------------------------------------
def permanova_phase(Y: pd.DataFrame, phases: list[str], permutations: int = 999) -> dict:
    """PERMANOVA on Bray-Curtis distance, testing Phase as the grouping factor."""
    D = beta_diversity('braycurtis', Y.values, ids=list(Y.index))
    grouping = pd.Series(phases, index=Y.index, name='Phase')
    res = skbio_permanova(D, grouping=grouping, permutations=permutations)
    # Convert pseudo-F to R² via SS_among / SS_total (PERMANOVA's variance partition)
    a = res['number of groups']
    n = res['sample size']
    f = res['test statistic']
    # R² from F: R² = F * (a-1) / (F * (a-1) + (n - a))
    r2 = float(f * (a - 1) / (f * (a - 1) + (n - a))) if (n - a) > 0 else float('nan')
    return {
        'F': float(f),
        'R2': r2,
        'p': float(res['p-value']),
        'n_perm': int(res['number of permutations']),
        'n_samples': int(n),
        'n_groups': int(a),
    }


def dbrda_perm_test(Y: pd.DataFrame, X: pd.DataFrame, *, n_perm: int = 999,
                    metric: str = 'braycurtis', seed: int = 42) -> dict:
    """Whole-model PERMANOVA-style test for db-RDA: shuffle Y rows and recompute the
    constrained-to-total inertia ratio (≅ vegan's anova(capscale, permutations=)).
    """
    obs = dbrda(Y, X, metric=metric)
    obs_ratio = obs['prop_constrained']
    rng = np.random.default_rng(seed)
    n = Y.shape[0]
    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(n)
        Y_perm = Y.iloc[idx].copy()
        Y_perm.index = Y.index
        null[i] = dbrda(Y_perm, X, metric=metric)['prop_constrained']
    p = (np.sum(null >= obs_ratio) + 1) / (n_perm + 1)
    return {
        'F_like_ratio': float(obs_ratio),
        'p': float(p),
        'n_perm': n_perm,
        'n_samples': n,
        'n_predictors': X.shape[1],
        'constrained_inertia': float(obs['constrained_inertia']),
        'total_inertia': float(obs['total_inertia']),
    }


# ----- biplot --------------------------------------------------------------
def plot_dbrda(res: dict, *, stages: list[str], sample_ids: list[str], title: str,
               arrow_color: str, out_path: Path, custom_arrow_labels: bool = False,
               ylim: tuple[float, float] | None = None,
               annotation: str | None = None, xs_mode: bool = False,
               arrow_label_offsets: dict[str, tuple[float, float]] | None = None) -> None:
    sites = res['sites']
    species = res['species'].copy()
    biplot = res['biplot'].copy()

    # Rescale species and biplot vectors to fit ~80% of the site half-range
    # (mirrors vegan's mul.arrow heuristic for scaling = 2 plots)
    site_half = max(
        sites['CAP1'].abs().max(),
        sites['CAP2'].abs().max(),
        1e-9,
    )
    sp_half = max(species[['CAP1', 'CAP2']].abs().max().max(), 1e-9)
    bp_half = max(biplot[['CAP1', 'CAP2']].abs().max().max(), 1e-9)
    sp_scale = 0.8 * site_half / sp_half
    bp_scale = 0.8 * site_half / bp_half
    species[['CAP1', 'CAP2']] *= sp_scale
    biplot[['CAP1', 'CAP2']] *= bp_scale

    fig, ax = plt.subplots(figsize=(10.7, 6.7))

    # 1. species labels — colored by GAO/PAO functional category (same scheme as heatmap axis labels)
    #    xs_mode renders every species marker as a literal "X" instead of its name
    for sp, row in species.iterrows():
        col = species_label_color(sp)
        fw = 'bold' if col != SPECIES_DEFAULT_COLOR else 'normal'
        ax.text(row['CAP1'], row['CAP2'], 'X' if xs_mode else sp, color=col,
                ha='center', va='center', fontsize=9 * 1.2 if xs_mode else 9, fontweight=fw)

    # 2. constraint arrows — layered above everything else, with bordered labels
    for var, row in biplot.iterrows():
        x, y = row['CAP1'], row['CAP2']
        ax.annotate('', xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2),
                    zorder=10)
        # default placement: out 10%
        x_lab, y_lab = x * 1.1, y * 1.1
        # match R's custom-label nudges for environment plot
        if custom_arrow_labels:
            if 'DO peak' in var or 'aerobic HRT' in var:
                x_lab, y_lab = x * 1.05, y - 0.07
            elif 'C/N' in var:
                x_lab, y_lab = x + 0.12, y
        # explicit per-variable label offsets (substring match) to separate
        # overlapping labels (e.g. peakN2O vs P removal in the performance panel)
        if arrow_label_offsets:
            for key, (dx, dy) in arrow_label_offsets.items():
                if key in var:
                    x_lab, y_lab = x_lab + dx, y_lab + dy
                    break
        ax.text(x_lab, y_lab, var, color=arrow_color, fontsize=12,
                ha='center', va='center', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                          edgecolor=arrow_color, linewidth=1.0, alpha=0.92),
                zorder=11)

    # 3. sample dots, colored by the dataset's own Phase assignment (stages is
    #    aligned to sample_ids). NOTE: the R script's B17–B60 stage_from_id map
    #    does NOT match this dataset (B32–B108) and mislabels Phase I/II, so we
    #    use stages directly.
    #    xs_mode: 20% larger circles (area ×1.44) with a thin black stroke
    cols = [STAGE_COLORS.get(s, '0.5') for s in stages]
    dot_size = 40 * 1.44 if xs_mode else 40
    ax.scatter(sites['CAP1'], sites['CAP2'], c=cols, s=dot_size, marker='o', alpha=0.9,
               edgecolors='black' if xs_mode else 'none',
               linewidths=0.6 if xs_mode else 0, zorder=3)

    # 4. axes, legend, title
    ax.axhline(0, color='lightgray', lw=0.8, zorder=0)
    ax.axvline(0, color='lightgray', lw=0.8, zorder=0)
    ax.set_xlabel(f'CAP1 ({100 * res["eig_constrained"].iloc[0] / res["total_inertia"]:.1f}%)')
    ax.set_ylabel(f'CAP2 ({100 * res["eig_constrained"].iloc[1] / res["total_inertia"]:.1f}%)' if len(res['eig_constrained']) > 1 else 'CAP2')
    ax.set_title(title)
    # set axes a bit beyond the union of site/species/biplot points so labels fit
    all_x = pd.concat([sites['CAP1'], species['CAP1'], biplot['CAP1']])
    all_y = pd.concat([sites['CAP2'], species['CAP2'], biplot['CAP2']])
    pad_x = 0.18 * (all_x.max() - all_x.min())
    pad_y = 0.18 * (all_y.max() - all_y.min())
    ax.set_xlim(all_x.min() - pad_x, all_x.max() + pad_x)
    ax.set_ylim(all_y.min() - pad_y, all_y.max() + pad_y)
    if ylim is not None:
        ax.set_ylim(*ylim)

    legend_handles = [
        plt.Line2D([0], [0], marker='s', color='none',
                   markerfacecolor=STAGE_COLORS[s], markeredgecolor='none',
                   markersize=10, label=f'Stage {s}')
        for s in ['I', 'II', 'III', 'IV', 'V']
    ]
    ax.legend(handles=legend_handles, loc='upper right', frameon=False)

    if annotation:
        ax.text(0.02, 0.02, annotation, transform=ax.transAxes,
                fontsize=9, va='bottom', ha='left', family='monospace',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          edgecolor='lightgray', alpha=0.92))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out_path}')

    # For the "X"-marker figures, export the genera each X represents, with the
    # functional category and the plotted (scaled) coordinates so each X on the
    # figure can be matched back to its genus.
    if xs_mode:
        sp_out = species[['CAP1', 'CAP2']].copy()
        sp_out.columns = ['plot_CAP1', 'plot_CAP2']
        sp_out.insert(0, 'functional_category', [species_category(g) for g in sp_out.index])
        sp_out.index.name = 'genus'
        sp_out = sp_out.sort_values(['functional_category', 'genus'])
        csv_path = out_path.with_name(out_path.stem + '_species.csv')
        sp_out.to_csv(csv_path)
        print(f'wrote {csv_path} ({len(sp_out)} genera)')


# ----- fallback: build stage-level Performance/Environmental from repo CSVs
PERF_KPIS = [
    'peakN2O [mg/L]',
    'specific denitrification rates [mg NO2–N g−1\xa0VSS−1 h−1]',
    'N removal (ppm) [N-ppn]',
    'P removal [P%]',
]
# Display-label overrides for the performance panel arrows (units stripped, etc.)
PERF_LABELS = {
    'specific denitrification rates [mg NO2–N g−1\xa0VSS−1 h−1]': 'specific denitrification rates',
}

# Curated "operational drivers" set replacing the old all-of-ValueEnviromental2 env
# panel. (source CSV, raw column name, display label) — paper 2's main influences.
INFLUENCE_VARS = [
    ('perf', 'COD:N',                'C/N'),
    ('perf', 'N:P',                  'N/P'),
    ('env',  'DO_Avg_1350_1410',     'DO avg'),
    ('perf', 'A_time [minutes]',     'cumulative aeration'),
    ('perf', 'Acetate ppm [mg/L]',   'Acetate'),
    ('perf', 'Propionate ppm [mg/L]','Propionate'),
    ('perf', 'N_Ax-1 [mg/L]',        'N_Ax-1'),
    ('perf', 'N_Ax-2 [mg/L]',        'N_Ax-2'),
]
PHASE_RANGES_FILE = ROOT / 'phase_day_ranges.json'
PERF_CSV = ROOT / 'performance_data_with_sample_ids.csv'
ENV_CSV  = ROOT / 'ValueEnviromental2_manual clean.csv'


def _abundance_day_to_phase(d: float, phase_ranges: dict[str, dict]) -> str | None:
    if pd.isna(d):
        return None
    for ph, span in phase_ranges.items():
        if span['start'] <= d <= span['end']:
            return ph
    return None


def build_perf_sheet_from_csv() -> pd.DataFrame:
    """Per-stage means of the 6 KPI columns, returned with a 'Stage' key."""
    phase_ranges = json.loads(PHASE_RANGES_FILE.read_text())
    perf = pd.read_csv(PERF_CSV, low_memory=False)
    cols = [c for c in PERF_KPIS if c in perf.columns]
    if not cols:
        raise RuntimeError(f'none of the expected KPIs {PERF_KPIS!r} are in {PERF_CSV}')
    perf['Phase'] = perf['abundance_day'].apply(lambda d: _abundance_day_to_phase(d, phase_ranges))
    sub = perf[perf['Phase'].isin(['I', 'II', 'III', 'IV', 'V'])][['Phase'] + cols].copy()
    out = sub.groupby('Phase').mean(numeric_only=True).reset_index().rename(columns={'Phase': 'Stage'})
    return out


def build_perf_X_per_sample(sample_ids: list[str], stages: list[str]) -> pd.DataFrame:
    """Per-sample Performance X.

    Performance KPIs in this dataset are sparse — many are recorded only at
    end-of-stage time points. Strategy: use each sample's directly-measured KPI
    value when present; otherwise impute with that sample's phase median so the
    samples still contribute to the dbRDA. Columns whose phase-medians are all
    NaN (no measurement in any phase) are dropped.
    """
    perf = pd.read_csv(PERF_CSV, low_memory=False)
    cols = [c for c in PERF_KPIS if c in perf.columns]
    if not cols:
        raise RuntimeError(f'none of the expected KPIs {PERF_KPIS!r} are in {PERF_CSV}')
    perf = perf[perf['sample_id'].isin(sample_ids)][['sample_id'] + cols].copy()
    perf = perf.set_index('sample_id').reindex(sample_ids)
    perf = perf.apply(pd.to_numeric, errors='coerce')

    phase_series = pd.Series(stages, index=sample_ids)
    phase_median = perf.groupby(phase_series).median()
    overall_median = perf.median()
    n_pre = perf.isna().sum().to_dict()
    for col in perf.columns:
        col_overall = overall_median[col]
        for sid, ph in zip(sample_ids, stages):
            if pd.isna(perf.at[sid, col]):
                v = phase_median.at[ph, col] if ph in phase_median.index else float('nan')
                if pd.isna(v):
                    v = col_overall  # fall back to global median if this phase had no measurements
                perf.at[sid, col] = v
    n_post = perf.isna().sum().to_dict()
    imputed = {c: n_pre[c] - n_post[c] for c in perf.columns if n_pre[c] > 0}
    if imputed:
        print(f'  perf: phase-median (overall-median fallback) imputation filled per-column: {imputed}')

    # Drop any column still NaN (would only happen if the KPI is never measured at all)
    bad_cols = [c for c in perf.columns if perf[c].isna().any()]
    if bad_cols:
        print(f'  perf: dropping columns that are never measured: {bad_cols}')
        perf = perf.drop(columns=bad_cols)
    return perf


def build_influence_X_per_sample(sample_ids: list[str], stages: list[str]) -> pd.DataFrame:
    """Per-sample operational-drivers matrix (paper-2 curated influences).

    Pulls the columns listed in INFLUENCE_VARS from their respective source CSV
    (performance vs ValueEnviromental2). Env columns are looked up on each
    sample's date via the 7-day-gap-interpolated daily env series. Missing
    values are imputed with phase median (overall median fallback) so all
    samples survive ``complete.cases``.
    """
    perf = pd.read_csv(PERF_CSV, low_memory=False)
    perf['date'] = pd.to_datetime(perf['date'], errors='coerce')
    perf = perf[perf['sample_id'].isin(sample_ids)]
    perf_indexed = perf.set_index('sample_id').reindex(sample_ids)
    sid_dates = perf.dropna(subset=['date']).drop_duplicates('sample_id')
    sample_to_date = dict(zip(sid_dates['sample_id'], sid_dates['date']))

    env = pd.read_csv(ENV_CSV)
    env['Date_key'] = pd.to_datetime(env['Date_key'], format='%m/%d/%Y', errors='coerce')
    env = env[env['Date_key'].notna()].set_index('Date_key').sort_index()
    drop_cols = [c for c in env.columns if c == 'Date' or c.startswith('Unnamed')]
    env = env.drop(columns=drop_cols, errors='ignore')
    env = env.apply(pd.to_numeric, errors='coerce').select_dtypes(include='number')
    env_daily = env.resample('D').mean().interpolate(method='time', limit=7, limit_area='inside')

    data: dict[str, pd.Series] = {}
    for source, src_col, display in INFLUENCE_VARS:
        if source == 'perf':
            if src_col not in perf_indexed.columns:
                print(f'  [warn] perf column missing: {src_col!r}')
                continue
            data[display] = pd.to_numeric(perf_indexed[src_col], errors='coerce')
        else:  # env
            if src_col not in env_daily.columns:
                print(f'  [warn] env column missing: {src_col!r}')
                continue
            data[display] = pd.Series(
                [env_daily.at[sample_to_date[s], src_col]
                 if (s in sample_to_date and sample_to_date[s] in env_daily.index) else np.nan
                 for s in sample_ids],
                index=sample_ids,
            )
    X = pd.DataFrame(data, index=sample_ids)

    # Impute: phase median → overall median
    phase_series = pd.Series(stages, index=sample_ids)
    phase_median = X.groupby(phase_series).median()
    overall_median = X.median()
    n_pre = X.isna().sum().to_dict()
    for col in X.columns:
        for sid, ph in zip(sample_ids, stages):
            if pd.isna(X.at[sid, col]):
                v = phase_median.at[ph, col] if ph in phase_median.index else float('nan')
                if pd.isna(v):
                    v = overall_median[col]
                X.at[sid, col] = v
    imputed = {c: n_pre[c] - X[c].isna().sum() for c in X.columns if n_pre[c] > 0}
    if imputed:
        print(f'  influences: phase-median (overall-median fallback) imputed: {imputed}')
    bad = [c for c in X.columns if X[c].isna().any()]
    if bad:
        print(f'  influences: dropping never-measured cols: {bad}')
        X = X.drop(columns=bad)
    return X


def build_env_X_per_sample(sample_ids: list[str]) -> pd.DataFrame:
    """Per-sample (not per-stage) Environmental X matrix.

    For each abundance sample, look up the ValueEnviromental2 row matching its
    sampling date (or interpolate within a 7-day gap, matching the operational
    correlation pipeline). Columns or samples with NaN are then trimmed.
    """
    env = pd.read_csv(ENV_CSV)
    env['Date_key'] = pd.to_datetime(env['Date_key'], format='%m/%d/%Y', errors='coerce')
    env = env[env['Date_key'].notna()].set_index('Date_key').sort_index()
    drop_cols = [c for c in env.columns if c == 'Date' or c.startswith('Unnamed')]
    env = env.drop(columns=drop_cols, errors='ignore')
    env = env.apply(pd.to_numeric, errors='coerce').select_dtypes(include='number')
    # daily index + 7-day gap-capped interpolation (matches correlation pipeline)
    env_daily = env.resample('D').mean().interpolate(method='time', limit=7, limit_area='inside')

    perf = pd.read_csv(PERF_CSV, low_memory=False)
    perf['date'] = pd.to_datetime(perf['date'], errors='coerce')
    sid_rows = perf[perf['sample_id'].isin(sample_ids) & perf['date'].notna()]
    sample_date = dict(zip(sid_rows['sample_id'], sid_rows['date']))

    rows = {}
    for sid in sample_ids:
        if sid in sample_date and sample_date[sid] in env_daily.index:
            rows[sid] = env_daily.loc[sample_date[sid]]
        else:
            rows[sid] = pd.Series(np.nan, index=env_daily.columns)
    X = pd.DataFrame(rows).T.reindex(sample_ids)
    # Drop columns that are entirely NaN (e.g., the string time-of-day cols)
    all_nan_cols = [c for c in X.columns if X[c].isna().all()]
    if all_nan_cols:
        print(f'  env: dropping all-NaN columns: {all_nan_cols}')
        X = X.drop(columns=all_nan_cols)
    # Downstream main() drops samples with any remaining NaN via complete.cases
    # (those are typically samples whose date falls outside the env data range)
    n_partial = X.isna().any(axis=1).sum()
    if n_partial:
        print(f'  env: {n_partial} samples have NaN in at least one env var '
              f'(typically outside env data range); will drop via complete.cases')
    return X


def build_env_sheet_from_csv() -> pd.DataFrame:
    """Per-stage means of every numeric env var in ValueEnviromental2, keyed by 'Stage'.

    Joins env dates to abundance_day via linear interpolation against the
    sample-day map in performance_data_with_sample_ids.csv, then averages within
    each phase's day window.
    """
    phase_ranges = json.loads(PHASE_RANGES_FILE.read_text())
    env = pd.read_csv(ENV_CSV)
    env['Date_key'] = pd.to_datetime(env['Date_key'], format='%m/%d/%Y', errors='coerce')
    env = env[env['Date_key'].notna()].set_index('Date_key').sort_index()
    drop_cols = [c for c in env.columns if c == 'Date' or c.startswith('Unnamed')]
    env = env.drop(columns=drop_cols, errors='ignore')
    env = env.apply(pd.to_numeric, errors='coerce').select_dtypes(include='number')

    perf = pd.read_csv(PERF_CSV, low_memory=False)
    perf['date'] = pd.to_datetime(perf['date'], errors='coerce')
    perf = perf.dropna(subset=['date', 'abundance_day']).sort_values('date')

    # Linear interpolation of abundance_day for every env date (within the perf range)
    perf_dates_ns = perf['date'].astype('int64').values
    perf_days = perf['abundance_day'].values
    env_dates_ns = env.index.astype('int64').values
    env['abundance_day'] = np.interp(
        env_dates_ns, perf_dates_ns, perf_days,
        left=np.nan, right=np.nan,
    )
    env['Phase'] = env['abundance_day'].apply(lambda d: _abundance_day_to_phase(d, phase_ranges))
    sub = env[env['Phase'].isin(['I', 'II', 'III', 'IV', 'V'])].drop(columns=['abundance_day'])
    out = sub.groupby('Phase').mean(numeric_only=True).reset_index().rename(columns={'Phase': 'Stage'})
    # Drop columns with any NaN at the stage level (e.g., the ANO_*_Time_C* string
    # columns that parsed as NaN, or env vars with no data in some phase)
    bad_cols = [c for c in out.columns if c != 'Stage' and out[c].isna().any()]
    if bad_cols:
        print(f'  dropping env cols with NaN at stage level: {bad_cols}')
        out = out.drop(columns=bad_cols)
    return out


# ----- broadcast a stage-level sheet to per-sample X ----------------------
def per_sample_X(stage_sheet: pd.DataFrame, sample_ids: list[str], stages: list[str]) -> pd.DataFrame:
    """Join a Stage-keyed predictor sheet to per-sample rows."""
    sheet = stage_sheet.copy()
    sheet['Stage'] = sheet['Stage'].astype(str)
    samp = pd.DataFrame({'sample': sample_ids, 'Stage': [str(s) for s in stages]})
    merged = samp.merge(sheet, on='Stage', how='left')
    X = merged.drop(columns=['sample', 'Stage']).apply(pd.to_numeric, errors='coerce')
    X.index = sample_ids
    return X


# ----- main ----------------------------------------------------------------
def main() -> None:
    print(f'Long abundance file: {LONG_FILE}')
    Y, sample_ids, stages = build_abundance_matrix(LONG_FILE)
    print(f'abundance matrix Y: {Y.shape[0]} samples × {Y.shape[1]} genera')

    print('using per-sample X (not stage-mean broadcast) — each sample contributes '
          'its own row from performance_data_with_sample_ids.csv / ValueEnviromental2')

    N_PERM = 999

    # ----- db-RDA: abundance ~ performance ----------------------------------
    X_perf = build_perf_X_per_sample(sample_ids, stages)
    X_perf = X_perf.rename(columns=PERF_LABELS)
    keep_p = X_perf.notna().all(axis=1)
    Y_p = Y.loc[keep_p]
    X_p = X_perf.loc[keep_p]
    stages_p = [s for s, k in zip(stages, keep_p) if k]
    stages_p = [s for s, k in zip(stages, keep_p) if k]
    print(f'\nperformance dbRDA: {Y_p.shape[0]} samples, {X_p.shape[1]} predictors')
    res_p = dbrda(Y_p, X_p)
    print(f'  constrained inertia: {res_p["constrained_inertia"]:.4f} '
          f'({100 * res_p["prop_constrained"]:.1f}% of total {res_p["total_inertia"]:.4f})')
    print(f'  running {N_PERM}-permutation PERMANOVAs...')
    pm_phase_p = permanova_phase(Y_p, stages_p, permutations=N_PERM)
    pm_full_p  = dbrda_perm_test(Y_p, X_p, n_perm=N_PERM)
    print(f'  PERMANOVA(Phase ~ Bray-Curtis): F={pm_phase_p["F"]:.2f}, '
          f'R²={pm_phase_p["R2"]:.3f}, p={pm_phase_p["p"]:.4f} '
          f'({pm_phase_p["n_groups"]} groups, {pm_phase_p["n_samples"]} samples)')
    print(f'  dbRDA whole-model perm test: constrained/total={pm_full_p["F_like_ratio"]:.3f}, '
          f'p={pm_full_p["p"]:.4f} ({pm_full_p["n_predictors"]} predictors)')
    ann_p = (
        f'PERMANOVA(Phase):  F={pm_phase_p["F"]:.2f}  R²={pm_phase_p["R2"]:.2f}  p={pm_phase_p["p"]:.3f}\n'
        f'db-RDA model:      constrained={100*pm_full_p["F_like_ratio"]:.1f}%  p={pm_full_p["p"]:.3f}\n'
        f'permutations={N_PERM}'
    )
    # peakN2O and P removal arrows point in nearly the same direction, so their
    # default labels overlap; pull them apart horizontally and drop them below
    # the sample-dot cluster.
    PERF_LABEL_OFFSETS = {
        'peakN2O':   (-0.06, -0.040),
        'P removal': (+0.060, -0.022),
    }
    plot_dbrda(
        res_p, stages=stages_p, sample_ids=Y_p.index.tolist(),
        title='db-RDA: Top 10 Genera ~ Performance (sample-level X)',
        arrow_color='crimson', out_path=OUT_DIR / 'dbRDA_performance.png',
        annotation=ann_p, arrow_label_offsets=PERF_LABEL_OFFSETS,
    )
    # _Xs variant: 20%-larger black-edged dots, species rendered as "X" (same model)
    plot_dbrda(
        res_p, stages=stages_p, sample_ids=Y_p.index.tolist(),
        title='db-RDA: Top 10 Genera ~ Performance (sample-level X)',
        arrow_color='crimson', out_path=OUT_DIR / 'dbRDA_performance_Xs.png',
        annotation=ann_p, xs_mode=True, arrow_label_offsets=PERF_LABEL_OFFSETS,
    )

    # ----- db-RDA: abundance ~ operational drivers (paper-2 curated set) ---
    X_inf = build_influence_X_per_sample(sample_ids, stages)
    keep_e = X_inf.notna().all(axis=1)
    Y_e = Y.loc[keep_e]
    X_e = X_inf.loc[keep_e]
    stages_e = [s for s, k in zip(stages, keep_e) if k]
    print(f'\noperational-drivers dbRDA: {Y_e.shape[0]} samples, {X_e.shape[1]} predictors')
    res_e = dbrda(Y_e, X_e)
    print(f'  constrained inertia: {res_e["constrained_inertia"]:.4f} '
          f'({100 * res_e["prop_constrained"]:.1f}% of total {res_e["total_inertia"]:.4f})')
    print(f'  running {N_PERM}-permutation PERMANOVAs...')
    pm_phase_e = permanova_phase(Y_e, stages_e, permutations=N_PERM)
    pm_full_e  = dbrda_perm_test(Y_e, X_e, n_perm=N_PERM)
    print(f'  PERMANOVA(Phase ~ Bray-Curtis): F={pm_phase_e["F"]:.2f}, '
          f'R²={pm_phase_e["R2"]:.3f}, p={pm_phase_e["p"]:.4f} '
          f'({pm_phase_e["n_groups"]} groups, {pm_phase_e["n_samples"]} samples)')
    print(f'  dbRDA whole-model perm test: constrained/total={pm_full_e["F_like_ratio"]:.3f}, '
          f'p={pm_full_e["p"]:.4f} ({pm_full_e["n_predictors"]} predictors)')
    ann_e = (
        f'PERMANOVA(Phase):  F={pm_phase_e["F"]:.2f}  R²={pm_phase_e["R2"]:.2f}  p={pm_phase_e["p"]:.3f}\n'
        f'db-RDA model:      constrained={100*pm_full_e["F_like_ratio"]:.1f}%  p={pm_full_e["p"]:.3f}\n'
        f'permutations={N_PERM}'
    )
    plot_dbrda(
        res_e, stages=stages_e, sample_ids=Y_e.index.tolist(),
        title='db-RDA: Top 10 Genera ~ Operational drivers (paper-2 set)',
        arrow_color='crimson', out_path=OUT_DIR / 'dbRDA_environment.png',
        annotation=ann_e,
    )

    # ----- _Xs variant: drop N_Ax-1, recompute, restyle (X markers, bigger edged dots) ---
    X_e_xs = X_e.drop(columns=[c for c in X_e.columns if c == 'N_Ax-1'])
    print(f'\noperational-drivers dbRDA (_Xs, N_Ax-1 removed): {Y_e.shape[0]} samples, {X_e_xs.shape[1]} predictors')
    res_e_xs = dbrda(Y_e, X_e_xs)
    print(f'  running {N_PERM}-permutation PERMANOVAs...')
    pm_phase_e_xs = permanova_phase(Y_e, stages_e, permutations=N_PERM)
    pm_full_e_xs  = dbrda_perm_test(Y_e, X_e_xs, n_perm=N_PERM)
    print(f'  dbRDA whole-model perm test: constrained/total={pm_full_e_xs["F_like_ratio"]:.3f}, '
          f'p={pm_full_e_xs["p"]:.4f} ({pm_full_e_xs["n_predictors"]} predictors)')
    ann_e_xs = (
        f'PERMANOVA(Phase):  F={pm_phase_e_xs["F"]:.2f}  R²={pm_phase_e_xs["R2"]:.2f}  p={pm_phase_e_xs["p"]:.3f}\n'
        f'db-RDA model:      constrained={100*pm_full_e_xs["F_like_ratio"]:.1f}%  p={pm_full_e_xs["p"]:.3f}\n'
        f'permutations={N_PERM}'
    )
    plot_dbrda(
        res_e_xs, stages=stages_e, sample_ids=Y_e.index.tolist(),
        title='db-RDA: Top 10 Genera ~ Operational drivers (paper-2 set, N_Ax-1 removed)',
        arrow_color='crimson', out_path=OUT_DIR / 'dbRDA_environment_Xs.png',
        annotation=ann_e_xs, xs_mode=True,
    )


if __name__ == '__main__':
    main()
