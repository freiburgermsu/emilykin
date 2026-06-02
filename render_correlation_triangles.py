"""Standalone renderer for the 'abundance correlation one-triangle' figures.

Mirrors render_diff_heatmaps.py / render_abundance_heatmaps.py: reads pre-computed
JSON/CSV artifacts in the working directory and re-renders the lower-triangle Spearman
correlation clustermaps WITHOUT re-running the upstream marimo pipeline. Run with the
uv-managed venv:

    ~/Documents/py_venv/bin/python render_correlation_triangles.py

Produces, each in a non-root and a root-aggregated flavour:
  - abundance_correlatons_one_triangle.png        full series, ASV-level labels
  - abundance_correlatons_one_triangle_root.png   full series, ASVs summed by root name
  - abundance_correlatons_one_triangle_phase{P}.png        per phase, ASV-level
  - abundance_correlatons_one_triangle_phase{P}_root.png   per phase, summed by root

Differences from the original inline figures:
  - Proteobacteria classes are shortened in the Phylum legend (Alphaproteobacteria
    -> Alpha, Gammaproteobacteria -> Gamma), matching the other figure suites.
  - Root flavours sum the relative abundances of every ASV sharing a phylogenetic
    root name (Thauera.34 + Thauera.7 + ... -> Thauera) before correlating.
  - The thick BLACK boxes marking H2-operationally-correlated members are omitted:
    their source data (model_inputs/, ASVset_correlations_*H2*.json) is no longer on
    disk. GAO/PAO colored boxes are unaffected.
"""

import json
import os
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from numpy import nan, ones_like, triu
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent
os.chdir(REPO)

TAXONOMIC_LEVELS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
GENUS_DEPTH = TAXONOMIC_LEVELS.index('Genus')
DEFAULT_COLOR = 'lightgray'
MIN_MAX_PCT = 0.005          # max relative abundance threshold (0.5%)
ZERO_LEVEL = 1e-05
EXCLUDE_DAYS_BELOW = 15      # per-phase: drop the reactor's first ~2 weeks
EXTRA_MAIN_ORG = 'Methanobacteriaceae.1'

GAOs_PAOs = json.load(open(REPO / 'gao_pao_categories.json'))
ITX = json.load(open(REPO / 'iterativeID_taxonomy.json'))
ICM = json.load(open(REPO / 'iterativeID_color_map.json'))
ILV = json.load(open(REPO / 'iterativeID_levels.json'))
PROTEO_CLASS_COLOR = json.load(open(REPO / 'proteo_class_color.json'))
PROTEO_BASE = json.load(open(REPO / 'phylum_base_overrides.json'))['Proteobacteria']

GENERA_COLOR_MAP = {ID.split('.')[0]: v for ID, v in ICM.items()}
INV_GAO_PAO = {v: k for k, vs in GAOs_PAOs.items() for v in vs}
LABEL_COLOR = {'GAOs': 'green', 'Putative GAOs': 'mediumseagreen',
               'PAOs': 'blue', 'Putative PAOs': 'cornflowerblue',
               'Other PHA storing potential+ function': 'red'}
BOX_COLOR = {'GAOs': 'green', 'Putative GAOs': 'mediumseagreen',
             'PAOs': 'blue', 'Putative PAOs': 'cornflowerblue'}

ID_LEVELS = {}
for _k, _v in ILV.items():
    ID_LEVELS[_k] = _v
    ID_LEVELS.setdefault(_k.split('.')[0], _v)

# root organism ID -> any constituent ASV's taxonomy dict
ROOT_TO_TAXA = {}
for _ID, _taxa in ITX.items():
    ROOT_TO_TAXA.setdefault(_ID.split('.')[0], _taxa)


def color_lookup(idx):
    if idx in ICM:
        return ICM[idx]
    if idx in GENERA_COLOR_MAP:
        return GENERA_COLOR_MAP[idx]
    return DEFAULT_COLOR


def taxonomy_string(org):
    taxa = ITX.get(org) or ROOT_TO_TAXA.get(org)
    if taxa is None:
        return f'Unknown|{org}'
    return '|'.join(str(taxa[l]) for l in TAXONOMIC_LEVELS[:GENUS_DEPTH + 1])


def spearman_matrix(abund):
    """Square Spearman correlation matrix over the columns of `abund`
    (rows = samples). NaN correlations (constant columns) become 0."""
    cols = list(abund.columns)
    n = len(cols)
    arr = abund.values.astype(float)
    M = np.ones((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            c, _ = spearmanr(arr[:, i], arr[:, j])
            if np.isnan(c):
                c = 0.0
            M[i, j] = M[j, i] = c
    return pd.DataFrame(M, index=cols, columns=cols)


def aggregate_to_root(abund):
    """Sum every ASV column sharing a phylogenetic root name into one column."""
    root = abund.copy()
    root.columns = [c.split('.')[0] for c in root.columns]
    return root.T.groupby(level=0).sum().T


def fdr_cooccurring(abund):
    """BH-FDR over Spearman p-values of all co-occurring pairs; returns the set
    of organisms that appear in at least one passing pair."""
    presence = (abund > ZERO_LEVEL).astype(int)
    cooc = defaultdict(int)
    for sample in presence.itertuples(index=False):
        present = [c for c, v in zip(presence.columns, sample) if v]
        for pair in combinations(sorted(present), 2):
            cooc[pair] += 1
    pdat = []
    for (a, b), _ct in cooc.items():
        r, p = spearmanr(abund[a], abund[b])
        if np.isnan(r):
            continue
        pdat.append((a, b, p))
    if not pdat:
        return set()
    pvals = np.array([t[2] for t in pdat])
    reject, _, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    orgs = set()
    for i, (a, b, _p) in enumerate(pdat):
        if reject[i]:
            orgs.add(a)
            orgs.add(b)
    return orgs


def _style_labels(ax):
    """Bold/color GAO-PAO organism tick labels and italicize genus-level labels."""
    for getter in (ax.get_yticklabels, ax.get_xticklabels):
        for lbl in getter():
            text = lbl.get_text()
            for org, cat in INV_GAO_PAO.items():
                if org in text:
                    lbl.set_fontweight('bold')
                    col = LABEL_COLOR.get(cat)
                    if col:
                        lbl.set_color(col)
            if ID_LEVELS.get(text) == 'Genus':
                lbl.set_fontstyle('italic')


def _phylum_legend_handles(corr_matrix, taxonomies):
    """Phylum legend handles with Proteobacteria split into shortened class names."""
    phylum_color = {}
    for idx in corr_matrix.index:
        parts = str(taxonomies.get(idx, '')).split('|')
        if len(parts) < 2:
            continue
        phylum = parts[1]
        if phylum in ('None', '', 'Unknown'):
            continue
        if phylum == 'Proteobacteria':
            phylum_color['Proteobacteria'] = PROTEO_BASE
            continue
        phylum_color.setdefault(phylum, color_lookup(idx))
    archaea_markers = ('archaeo', 'halobacterota', 'methanobacteriota')

    def is_archaea(p):
        return any(m in p.lower() for m in archaea_markers)

    archaea = sorted(p for p in phylum_color if is_archaea(p))
    bacteria = sorted(p for p in phylum_color if not is_archaea(p))
    handles = []
    if archaea:
        handles.append(mpatches.Patch(color='none', label=r'$\bf{Archaea}$'))
        handles += [mpatches.Patch(facecolor=phylum_color[p], label=p) for p in archaea]
    if bacteria:
        handles.append(mpatches.Patch(color='none', label=r'$\bf{Bacteria}$'))
        for p in bacteria:
            if p == 'Proteobacteria' and PROTEO_CLASS_COLOR:
                handles.append(mpatches.Patch(color='none', label=p))
                for cls, color in PROTEO_CLASS_COLOR.items():
                    cls_short = cls.replace('proteobacteria', '').replace('Proteobacteria', '')
                    handles.append(mpatches.Patch(facecolor=color, label=f'      {cls_short}'))
            else:
                handles.append(mpatches.Patch(facecolor=phylum_color[p], label=p))
    return handles


def render_one_triangle(corr_matrix, taxonomies, out_path):
    """Lower-triangle Spearman correlation clustermap with phylum color strips,
    GAO/PAO colored boxes, a shortened-Proteobacteria Phylum legend, and the
    colorbar tucked into the empty upper-right of the triangle.
    """
    if corr_matrix.shape[0] < 2:
        print(f'[{out_path}] not enough organisms ({corr_matrix.shape[0]}); skipped')
        return
    row_colors = pd.Series({i: color_lookup(i) for i in corr_matrix.index}, name='Phylum')
    col_colors = pd.Series({c: color_lookup(c) for c in corr_matrix.columns}, name='Phylum')

    cm = sns.clustermap(
        corr_matrix, row_colors=row_colors, col_colors=col_colors,
        cmap='coolwarm_r', center=0, figsize=(60, 70), dendrogram_ratio=(0.1, 0.2),
    )
    cm.figure.subplots_adjust(bottom=0.15, top=0.95)
    cm.ax_row_dendrogram.set_visible(False)
    cm.ax_col_dendrogram.set_visible(False)
    cm.ax_heatmap.yaxis.set_ticks_position('left')
    cm.ax_heatmap.yaxis.set_label_position('left')

    hm_pos = cm.ax_heatmap.get_position()
    fig_w = cm.figure.get_figwidth()
    fig_h = cm.figure.get_figheight()
    strip_w = strip_h = 0.015
    cm.ax_row_colors.set_position([hm_pos.x0 - strip_w, hm_pos.y0, strip_w, hm_pos.height])
    cm.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0 - strip_h, hm_pos.width, strip_h])
    cm.ax_heatmap.tick_params(axis='y', pad=strip_w * fig_w * 72 + 15)
    cm.ax_heatmap.tick_params(axis='x', pad=strip_h * fig_h * 72 + 15)

    cbar = cm.ax_cbar
    ticks = list(cbar.get_yticks())
    ticks[-1] = corr_matrix.max().max()
    ticks[0] = round(corr_matrix.min().min(), 2)
    cbar.set_yticks(ticks)
    cbar.set_yticklabels([f'{t:g}' for t in ticks], fontsize=45)
    cbar.set_xlabel('Spearman Correlation', fontsize=45, labelpad=30)

    cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(), fontsize=40, rotation=0)
    cm.ax_heatmap.set_xticklabels(
        cm.ax_heatmap.get_xticklabels(), fontsize=40, rotation=60, ha='right', rotation_mode='anchor')
    _style_labels(cm.ax_heatmap)

    # Mask the upper triangle so only the lower triangle is drawn.
    d_row = cm.dendrogram_row.reordered_ind
    d_col = cm.dendrogram_col.reordered_ind
    df_re = corr_matrix.iloc[d_row, d_col]
    mask = triu(ones_like(df_re, dtype=bool), k=1)
    mesh = cm.ax_heatmap.collections[0]
    arr = mesh.get_array().reshape(df_re.shape)
    arr[mask] = nan
    mesh.set_array(arr.ravel())

    heatmap_pos = cm.ax_heatmap.get_position()
    cm.ax_cbar.set_position([heatmap_pos.x1 - 0.2, heatmap_pos.y0 + 0.4, 0.06, heatmap_pos.height / 5])

    # GAO/PAO colored boxes (one row strip + one column strip per organism).
    n_corr = len(corr_matrix.index)
    for _id in corr_matrix.index:
        cat = INV_GAO_PAO.get(_id.split('.')[0])
        if cat not in BOX_COLOR:
            continue
        color = BOX_COLOR[cat]
        rp = d_row.index(corr_matrix.index.get_loc(_id))
        cp = d_col.index(corr_matrix.columns.get_loc(_id))
        cm.ax_heatmap.add_patch(mpatches.Rectangle(
            (0, rp), rp + 1, 1, linewidth=4, edgecolor=color, facecolor='none', clip_on=False))
        cm.ax_heatmap.add_patch(mpatches.Rectangle(
            (cp, n_corr), 1, -(n_corr - cp), linewidth=4, edgecolor=color, facecolor='none', clip_on=False))

    handles = _phylum_legend_handles(corr_matrix, taxonomies)
    cm.figure.legend(
        handles=handles, title='Phylum', title_fontsize=50, fontsize=40,
        loc='upper right', bbox_to_anchor=(0.65, 0.7), frameon=True,
        borderaxespad=0.5, handlelength=1.5, handletextpad=0.6)
    cm.ax_heatmap.set_xlabel('Member ASVs', fontsize=50)
    cm.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)

    cm.figure.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(cm.figure)
    print(f'wrote {out_path}  ({corr_matrix.shape[0]} organisms)')


def build_taxonomies(orgs):
    return {o: taxonomy_string(o) for o in orgs}


def main():
    sample_days = json.load(open('sample_days.json'))
    phase_day_ranges = json.load(open('phase_day_ranges.json'))

    abundances = pd.read_csv('abundances.csv', header=0).set_index('sample')

    # ---- Main full-series figure ----
    sig = {str(x) for x in np.load('FDR_passing_pairs.npy', allow_pickle=True)}
    sig.add(EXTRA_MAIN_ORG)

    # non-root: significant ASVs with max relative abundance >= 0.5%
    asv_cols = [c for c in abundances.columns if c in sig and abundances[c].max() >= MIN_MAX_PCT]
    corr_asv = spearman_matrix(abundances[asv_cols])
    render_one_triangle(corr_asv, build_taxonomies(asv_cols),
                        'abundance_correlatons_one_triangle.png')

    # root: sum all ASVs per root, keep roots with a significant member and max >= 0.5%
    sig_roots = {c.split('.')[0] for c in sig}
    root_abund = aggregate_to_root(abundances)
    root_cols = [r for r in root_abund.columns if r in sig_roots and root_abund[r].max() >= MIN_MAX_PCT]
    corr_root = spearman_matrix(root_abund[root_cols])
    render_one_triangle(corr_root, build_taxonomies(root_cols),
                        'abundance_correlatons_one_triangle_root.png')

    # ---- Per-phase figures ----
    for phase, span in phase_day_ranges.items():
        lo = max(span['start'], EXCLUDE_DAYS_BELOW)
        keep_samples = [s for s in abundances.index if lo <= int(sample_days[s]) <= span['end']]
        abp = abundances.loc[keep_samples]
        abp = abp.loc[:, (abp.fillna(0) > 0).any(axis=0)]
        if abp.shape[0] < 3 or abp.shape[1] < 2:
            print(f'[phase {phase}] too few samples/organisms; skipped')
            continue

        # non-root
        po = fdr_cooccurring(abp)
        if po:
            sig_phase = [c for c in abp.columns if c in po]
            abc = abp[sig_phase].loc[:, abp[sig_phase].max() >= MIN_MAX_PCT]
            if abc.shape[1] >= 2:
                corr = spearman_matrix(abc)
                render_one_triangle(corr, build_taxonomies(abc.columns),
                                    f'abundance_correlatons_one_triangle_phase{phase}.png')

        # root: sum the phase table to root, then the same FDR + max filter
        abp_root = aggregate_to_root(abp)
        abp_root = abp_root.loc[:, (abp_root > 0).any(axis=0)]
        po_r = fdr_cooccurring(abp_root)
        if po_r:
            sig_r = [c for c in abp_root.columns if c in po_r]
            abc_r = abp_root[sig_r].loc[:, abp_root[sig_r].max() >= MIN_MAX_PCT]
            if abc_r.shape[1] >= 2:
                corr_r = spearman_matrix(abc_r)
                render_one_triangle(corr_r, build_taxonomies(abc_r.columns),
                                    f'abundance_correlatons_one_triangle_phase{phase}_root.png')


if __name__ == '__main__':
    main()
