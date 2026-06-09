"""Standalone renderer for the 'Top 10 ASVs (% abundance)' relative-abundance heatmaps.

Mirrors render_diff_heatmaps.py: reads the pre-computed JSON/CSV artifacts in the
working directory and re-renders the figure logic that otherwise lives inline in the
data_processing.py marimo cell, WITHOUT re-running the upstream pipeline. Run with the
uv-managed venv:

    ~/Documents/py_venv/bin/python render_abundance_heatmaps.py

Produces (each in a non-root and a root-aggregated flavour):
  - top_10_asvs_(%_abundance).png       top 10 ASVs/day, ASV-level labels (Thauera.34)
  - top_10_asvs_(%_abundance)_root.png   top 10 root organisms/day, abundances of all
                                         ASVs sharing a phylogenetic root name summed
                                         (Thauera.34 + Thauera.7 + ... -> Thauera)

Both flavours shorten the Proteobacteria classes in the Phylum legend
(Alphaproteobacteria -> Alpha, Gammaproteobacteria -> Gamma).
"""

import json
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
from numpy import inf, isnan, log10, nan

REPO = Path(__file__).resolve().parent
# Scripts now live in scripts/; the data artifacts sit in the repo root.
if not (REPO / 'abundances.csv').exists() and (REPO.parent / 'abundances.csv').exists():
    REPO = REPO.parent
os.chdir(REPO)
OUT_DIR = REPO / 'rel_ab_heatmaps'
OUT_DIR.mkdir(exist_ok=True)

TAXONOMIC_LEVELS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
GENUS_DEPTH = TAXONOMIC_LEVELS.index('Genus')
DEFAULT_COLOR = 'lightgray'
ZERO_LEVEL = 1e-05
TOP_NUM = 10

GAOs_PAOs = json.load(open(REPO / 'gao_pao_categories.json'))


def taxonomy_linkage(taxonomy_series):
    """Build a linkage matrix that exactly follows taxonomy hierarchy."""
    parsed = (taxonomy_series.fillna('unknown').astype(str)
              .str.split('[;|,]', regex=True)
              .apply(lambda lst: [x.strip() for x in lst or [] if x.strip()]))
    ranks_n = max((len(p) for p in parsed))
    if ranks_n == 0:
        ranks_n = 1
    parsed = parsed.apply(lambda lst: (lst + [''] * ranks_n)[:ranks_n])
    n = len(parsed)
    if n <= 1:
        return np.empty((0, 4), dtype=float)
    idx_to_pos = {idx: pos for pos, idx in enumerate(parsed.index)}
    linkage_rows = []
    next_cluster_id = n

    def build_subtree(indices, depth):
        nonlocal next_cluster_id
        if len(indices) == 1:
            return (idx_to_pos[indices[0]], 1)
        if depth >= ranks_n:
            merge_distance = 0.5
            cluster_id, total_count = (idx_to_pos[indices[0]], 1)
            for idx in indices[1:]:
                new_id = next_cluster_id
                next_cluster_id += 1
                linkage_rows.append([cluster_id, idx_to_pos[idx], merge_distance, total_count + 1])
                cluster_id = new_id
                total_count += 1
            return (cluster_id, total_count)
        groups = defaultdict(list)
        for idx in indices:
            rank_value = parsed.loc[idx][depth]
            if rank_value == '':
                rank_value = f'__unclassified_{idx}'
            groups[rank_value].append(idx)
        subclusters = []
        for rank_value, group_indices in groups.items():
            sub_id, sub_count = build_subtree(group_indices, depth + 1)
            subclusters.append((sub_id, sub_count))
        if len(subclusters) == 1:
            return subclusters[0]
        merge_distance = float(ranks_n - depth)
        cluster_id, total_count = subclusters[0]
        for sub_id, sub_count in subclusters[1:]:
            new_id = next_cluster_id
            next_cluster_id += 1
            linkage_rows.append([cluster_id, sub_id, merge_distance, total_count + sub_count])
            cluster_id = new_id
            total_count += sub_count
        return (cluster_id, total_count)

    build_subtree(list(parsed.index), depth=0)
    Z = np.array(linkage_rows, dtype=float)
    if len(Z) > 0:
        for i in range(1, len(Z)):
            if Z[i, 2] < Z[i - 1, 2]:
                Z[i, 2] = Z[i - 1, 2]
    return Z


def create_heatmap(df, taxonomies, title, *, suffix, simple=False,
                   iterativeID_color_map, genera_color_map,
                   proteo_class_color, proteo_base, iterativeID_levels,
                   phase_day_ranges):
    """Render the log10 relative-abundance heatmap: phylum row colors on the right,
    organism labels on the left, taxonomy dendrogram on the far right, vertical
    dashed phase delimiters, and a Phylum legend in the upper right.
    """
    new_cmap = LinearSegmentedColormap.from_list(
        'NewMap', [(0.0, 'aliceblue'), (0.25, 'lightblue'), (1.0, 'navy')])
    new_cmap.set_bad('aliceblue')
    vmin = df.min().min()
    vmax = df.max().max()
    vcenter = log10(0.1)
    vcenter = min(max(vcenter, vmin + 1e-06), vmax - 1e-06)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    def lookup(idx):
        if idx in iterativeID_color_map:
            return iterativeID_color_map[idx]
        if idx in genera_color_map:
            return genera_color_map[idx]
        return DEFAULT_COLOR
    row_colors = pd.Series({idx: lookup(idx) for idx in df.index}, name='Phylum')

    figsize = (50, 20)
    dendro_ratio = (0.2, 0.15)
    cm = sns.clustermap(
        df, row_colors=row_colors, cmap=new_cmap, norm=norm, clip_on=True,
        col_cluster=False, row_cluster=True, figsize=figsize,
        row_linkage=taxonomy_linkage(taxonomies), dendrogram_ratio=dendro_ratio,
    )
    for tick in cm.ax_row_colors.get_xticklabels():
        tick.set_fontsize(22)

    cbar = cm.ax_cbar
    log_ticks = np.array([vmin, log10(0.01), vcenter, log10(0.2), vmax])
    cbar.set_yticks(log_ticks)
    original_ticks = 10 ** log_ticks
    original_labels = [f'{x * 100:.1e}'.replace('e-0', 'E-') if x * 100 < 0.1 else round(x * 100)
                       for x in original_ticks]
    original_labels[0] = '0'
    cbar.set_yticklabels(original_labels, fontsize=22)
    cbar.set_yticks(log_ticks)
    cbar.set_xlabel('Rel. Abundance %', fontsize=18, labelpad=10)

    cm.ax_col_dendrogram.set_visible(False)
    cm.figure.subplots_adjust(bottom=0.15, top=0.95)
    cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(), fontsize=24, rotation=0)
    cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(), fontsize=24, rotation=70, ha='right')
    cm.ax_heatmap.set_xlabel('Days of operation', fontsize=32, labelpad=20)
    cm.ax_row_dendrogram.xaxis.set_visible(False)
    cm.ax_row_dendrogram.text(0.65, -0.02, 'Taxonomical tree', fontsize=24, ha='center',
                              transform=cm.ax_row_dendrogram.transAxes)

    # Layout: [row labels] [heatmap] [phylum colors] [dendrogram], cbar on the far left.
    cm.ax_heatmap.yaxis.tick_left()
    cm.ax_heatmap.yaxis.set_label_position('left')
    label_right = 0.22
    heatmap_right = 0.68
    dendro_left = 0.71
    dendro_right = 0.85
    hm_pos = cm.ax_heatmap.get_position()
    dend_pos = cm.ax_row_dendrogram.get_position()
    cm.ax_heatmap.set_position([label_right, hm_pos.y0, heatmap_right - label_right, hm_pos.height])
    cm.ax_row_dendrogram.set_position([dendro_left, dend_pos.y0, dendro_right - dendro_left, dend_pos.height])
    cm.ax_row_dendrogram.invert_xaxis()
    gap = 0.03
    col_dend_pos = cm.ax_col_dendrogram.get_position()
    cm.ax_col_dendrogram.set_position([label_right, col_dend_pos.y0 + gap, heatmap_right - label_right, col_dend_pos.height])
    top_pos = cm.ax_col_dendrogram.get_position()
    cbar_w = 0.05
    cbar_x = top_pos.x0 - cbar_w - 0.02
    cm.ax_cbar.set_position([cbar_x, top_pos.y0, cbar_w, top_pos.height])

    if simple:
        # Simplified labelling: only Ca_Accumulibacter (blue) and Ca_Competibacter
        # (green) are highlighted; every other label stays plain black.
        for label in cm.ax_heatmap.get_yticklabels():
            text = label.get_text()
            if 'Ca_Accumulibacter' in text:
                label.set_color('blue')
                label.set_fontweight('bold')
            elif 'Ca_Competibacter' in text:
                label.set_color('green')
                label.set_fontweight('bold')
    else:
        # GAO/PAO label coloring and italicization at the genus level.
        id_levels = {}
        for k, v in iterativeID_levels.items():
            id_levels[k] = v
            id_levels.setdefault(k.split('.')[0], v)
        inverted_gao_pao = {v: k for k, vs in GAOs_PAOs.items() for v in vs}
        for label in cm.ax_heatmap.get_yticklabels():
            text = label.get_text()
            for org, val in inverted_gao_pao.items():
                if org not in text:
                    continue
                label.set_fontweight('bold')
                if 'GAOs' in val:
                    label.set_color('green')
                    if 'Putative' in val:
                        label.set_color('mediumseagreen')
                elif 'PAOs' in val:
                    label.set_color('blue')
                    if 'Putative' in val:
                        label.set_color('cornflowerblue')
                elif 'PHA' in val:
                    label.set_color('red')
            if id_levels.get(text) == 'Genus':
                label.set_fontstyle('italic')

    # Phylum color strip to the RIGHT of the heatmap.
    hm_pos = cm.ax_heatmap.get_position()
    strip_w = 0.015
    cm.ax_row_colors.set_position([hm_pos.x1 + 0.005, hm_pos.y0, strip_w, hm_pos.height])
    cm.ax_heatmap.tick_params(axis='y', pad=4)

    # Build phylum legend (one row per phylum, plus Proteobacteria classes).
    phylum_color = {}
    for idx in df.index:
        parts = str(taxonomies.get(idx, '')).split('|')
        if len(parts) < 2:
            continue
        phylum = parts[1]
        if phylum in ('None', '', 'Unknown'):
            continue
        if phylum == 'Proteobacteria':
            phylum_color['Proteobacteria'] = proteo_base
            continue
        color = row_colors.get(idx)
        if color is None:
            continue
        phylum_color.setdefault(phylum, color)
    archaea_markers = ('archaeo', 'halobacterota', 'methanobacteriota')

    def is_archaea(p):
        return any(m in p.lower() for m in archaea_markers)

    archaea = sorted(p for p in phylum_color if is_archaea(p))
    bacteria = sorted(p for p in phylum_color if not is_archaea(p))
    handles = []
    if archaea:
        handles.append(Patch(color='none', label=r'$\bf{Archaea}$'))
        handles += [Patch(facecolor=phylum_color[p], label=p) for p in archaea]
    if bacteria:
        handles.append(Patch(color='none', label=r'$\bf{Bacteria}$'))
        for p in bacteria:
            if p == 'Proteobacteria' and proteo_class_color:
                handles.append(Patch(color='none', label=p))
                for cls, color in proteo_class_color.items():
                    cls_short = cls.replace('proteobacteria', '').replace('Proteobacteria', '')
                    handles.append(Patch(facecolor=color, label=f'      {cls_short}'))
            else:
                handles.append(Patch(facecolor=phylum_color[p], label=p))
    cm.figure.legend(
        handles=handles, title='Phylum', title_fontsize=20, fontsize=18,
        loc='upper right', bbox_to_anchor=(dendro_right, 1.0), frameon=True,
        borderaxespad=0.5, handlelength=1.5, handletextpad=0.6,
    )
    for spine in cm.ax_heatmap.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(1)

    # Vertical dashed phase delimiters at the lower bound of each phase (skip first).
    if phase_day_ranges:
        days_int = []
        for c in df.columns:
            try:
                days_int.append(int(c))
            except (ValueError, TypeError):
                days_int = None
                break
        if days_int:
            n_cols = len(days_int)

            def day_to_heatmap_x(d):
                if d <= days_int[0]:
                    return 0.5
                if d >= days_int[-1]:
                    return n_cols - 0.5
                for i in range(n_cols - 1):
                    if days_int[i] <= d <= days_int[i + 1]:
                        if days_int[i + 1] == days_int[i]:
                            return i + 0.5
                        frac = (d - days_int[i]) / (days_int[i + 1] - days_int[i])
                        return i + 0.5 + frac
                return None

            fig = cm.figure
            hm_pos = cm.ax_heatmap.get_position()
            top_pos = cm.ax_col_dendrogram.get_position()
            y_bot = hm_pos.y0
            y_top = top_pos.y1
            fig.canvas.draw()
            lower_bounds = sorted({span['start'] for span in phase_day_ranges.values()})[1:]
            for day in lower_bounds:
                x_data = day_to_heatmap_x(day)
                if x_data is None:
                    continue
                x_disp = cm.ax_heatmap.transData.transform((x_data, 0))[0]
                x_fig = fig.transFigure.inverted().transform((x_disp, 0))[0]
                fig.add_artist(mlines.Line2D(
                    [x_fig, x_fig], [y_bot, y_top], transform=fig.transFigure,
                    color='black', linewidth=1.2, linestyle='--', alpha=0.6))

    out_path = OUT_DIR / f"{title.lower().replace(' ', '_')}{suffix}.png"
    cm.figure.savefig(out_path, bbox_inches='tight', dpi=300)
    print(f'wrote {out_path}')


def taxonomy_string(taxa):
    """'|'-joined taxonomy from Kingdom..Genus (matches the inline pipeline)."""
    return '|'.join(str(taxa[l]) for l in TAXONOMIC_LEVELS[:GENUS_DEPTH + 1])


def top10_frame(abundances, org_taxonomy):
    """Replicate the inline top-10-per-day selection.

    abundances: rows = days, cols = organism IDs, values = relative abundance.
    org_taxonomy: org ID -> '|'-joined taxonomy string.
    Returns (df, taxonomy_series) where df is log10(abundance) over the union of
    each day's top-10 organisms (rows = orgs, cols = days).
    """
    dic = {}
    for day, row in abundances.iterrows():
        for org, ab in row.items():
            if ab <= ZERO_LEVEL:
                continue
            dic.setdefault(day, {})
            dic[day][org] = dic[day].get(org, 0) + ab
    nonzero_per_day = {
        day: dict(sorted(org_dict.items(), key=lambda kv: kv[1], reverse=True))
        for day, org_dict in dic.items()
    }
    top10_all_days = []
    for org_dict in nonzero_per_day.values():
        top10_all_days.extend(list(org_dict.keys())[:TOP_NUM])
    top10_all_days = set(top10_all_days)
    top_per_day = {
        day: {org: log10(v) for org, v in org_dict.items() if org in top10_all_days}
        for day, org_dict in nonzero_per_day.items()
    }
    df = pd.DataFrame(top_per_day).astype(float).replace([inf, -inf], nan)
    taxonomy_series = pd.Series({idx: org_taxonomy.get(idx, f'Unknown|{idx}') for idx in df.index})
    return df, taxonomy_series


def main(simple_values=(False, True)):
    sample_days = json.load(open('sample_days.json'))
    iterativeID_color_map = json.load(open('iterativeID_color_map.json'))
    proteo_class_color = json.load(open('proteo_class_color.json'))
    proteo_base = json.load(open('phylum_base_overrides.json'))['Proteobacteria']
    iterativeID_taxonomy = json.load(open('iterativeID_taxonomy.json'))
    iterativeID_levels = json.load(open('iterativeID_levels.json'))
    phase_day_ranges = json.load(open('phase_day_ranges.json'))

    genera_color_map = {ID.split('.')[0]: v for ID, v in iterativeID_color_map.items()}

    abundances = pd.read_csv('abundances.csv').set_index('sample')
    abundances.index = [sample_days[col] for col in abundances.index]
    abundances = abundances.loc[sorted(abundances.index, key=int)]

    shared_kwargs = dict(
        iterativeID_color_map=iterativeID_color_map,
        genera_color_map=genera_color_map,
        proteo_class_color=proteo_class_color,
        proteo_base=proteo_base,
        iterativeID_levels=iterativeID_levels,
        phase_day_ranges=phase_day_ranges,
    )

    # ---- Non-root: top 10 ASVs/day, ASV-level labels ----
    asv_taxonomy = {org: taxonomy_string(iterativeID_taxonomy[org]) for org in abundances.columns}
    df_asv, tax_asv = top10_frame(abundances, asv_taxonomy)

    # ---- Root: sum all ASVs sharing a phylogenetic root name, then top 10/day ----
    root_abundances = abundances.copy()
    root_abundances.columns = [c.split('.')[0] for c in root_abundances.columns]
    root_abundances = root_abundances.T.groupby(level=0).sum().T
    root_taxonomy = {}
    for ID, taxa in iterativeID_taxonomy.items():
        root = ID.split('.')[0]
        if root in root_abundances.columns and root not in root_taxonomy:
            root_taxonomy[root] = taxonomy_string(taxa)
    df_root, tax_root = top10_frame(root_abundances, root_taxonomy)

    # `simple` flavours suffix '_simple'; the GAO/PAO-rich default has no suffix.
    for simple in simple_values:
        s = '_simple' if simple else ''
        create_heatmap(df_asv, tax_asv, 'Top 10 ASVs (% abundance)', suffix=s,
                       simple=simple, **shared_kwargs)
        create_heatmap(df_root, tax_root, 'Top 10 ASVs (% abundance)', suffix='_root' + s,
                       simple=simple, **shared_kwargs)


if __name__ == '__main__':
    main()
