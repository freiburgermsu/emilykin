"""Standalone renderer for the 'taxa above 1% max abundance' fold-change heatmaps.

Reads pre-computed JSON/CSV artifacts in the working directory; does NOT
regenerate any of the upstream pipeline outputs (color maps, taxonomy JSONs,
nonzero_per_day.json, etc.). Run with the codiffusion_bioreactor venv:

    /home/andrew/repos/codiffusion_bioreactor/.venv/bin/python render_diff_heatmaps.py

Produces:
  - taxa_above_1%_max_abundance_diff.png            log2 FC between phase pairs
  - taxa_above_1%_max_abundance_diff_innoculum.png  log2 FC per day vs the first recorded day
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
from scipy.cluster.hierarchy import linkage

REPO = Path(__file__).resolve().parent
os.chdir(REPO)

TAXONOMIC_LEVELS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
DEFAULT_COLOR = 'lightgray'
MIN_MAX_PCT = 0.01
PSEUDO = 1e-6
GENUS_DEPTH = TAXONOMIC_LEVELS.index('Genus')

GAOs_PAOs = {
    'PAOs': ['Ca_Accumulibacter', 'Tetrasphaera', 'Dechloromonas', 'Microlunatus', 'Azonexus', 'Ca_Phosphoribacter'],
    'GAOs': ['Ca_Competibacter', 'Defluviicoccus', 'Propionivibrio', 'Ca_Contendobacter'],
    'Putative PAOs': ['Ca_Obscuribacter', 'Thauera', 'Zoogloea', 'Paracoccus'],
    'Putative GAOs': ['Micropruina', 'Amaricoccus', 'Ca_Glycocaulis', 'Thauera'],
    'Other PHA storing potential+ function': ['Pseudomonas', 'Bacillus', 'Acinetobacter', 'Rhodocyclaceae'],
}


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


def _draw_phase_delimiters(clustermap, df, phase_day_ranges):
    """Vertical dashed lines at the lower bound of each phase (skipping the first)."""
    days_int = []
    for c in df.columns:
        try:
            days_int.append(int(c))
        except (ValueError, TypeError):
            return
    n_cols = len(days_int)

    def day_to_x(d):
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

    fig = clustermap.figure
    hm_pos = clustermap.ax_heatmap.get_position()
    y_bot = hm_pos.y0
    y_top = hm_pos.y1  # stop at the top edge of the heatmap; don't extend above it
    fig.canvas.draw()
    lower_bounds = sorted({span['start'] for span in phase_day_ranges.values()})[1:]
    for day in lower_bounds:
        x_data = day_to_x(day)
        if x_data is None:
            continue
        x_disp = clustermap.ax_heatmap.transData.transform((x_data, 0))[0]
        x_fig = fig.transFigure.inverted().transform((x_disp, 0))[0]
        fig.add_artist(mlines.Line2D(
            [x_fig, x_fig], [y_bot, y_top], transform=fig.transFigure,
            color='black', linewidth=1.2, linestyle='--', alpha=0.6))


def create_heatmap(df, taxonomies, title, *, mode, suffix,
                   iterativeID_color_map, genera_color_map,
                   proteo_class_color, proteo_base, iterativeID_levels,
                   phase_day_ranges):
    """Render a log2-fold-change heatmap with phylum row colors on the left,
    organism labels between heatmap and dendrogram, and a horizontal phylum
    legend centered at the top.

    mode='log2fc'      → columns are phase pairs (compact figure, no day delims)
    mode='log2fc_days' → columns are days (wide figure, vertical phase delimiters)
    """
    cmap = matplotlib.colormaps['coolwarm_r'].copy()
    cmap.set_bad('lightgray')
    abs_max = float(np.nanmax(np.abs(df.values))) if df.size else 1.0
    if not np.isfinite(abs_max) or abs_max == 0:
        abs_max = 1.0
    vmin, vmax = -abs_max, abs_max
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    def lookup(idx):
        if idx in iterativeID_color_map:
            return iterativeID_color_map[idx]
        if idx in genera_color_map:
            return genera_color_map[idx]
        return DEFAULT_COLOR
    row_colors = pd.Series({idx: lookup(idx) for idx in df.index}, name='Phylum')

    if mode == 'log2fc_days':
        figsize = (50, 20)
        dendro_ratio = (0.06, 0.15)
        heatmap_right = 0.72
        dendro_w = 0.04
        xlabel = 'Days of operation'
    else:
        figsize = (max(14, 1.5 * df.shape[1] + 8), 20)
        dendro_ratio = (0.08, 0.15)
        heatmap_right = 0.58
        dendro_w = 0.06
        xlabel = 'Phase comparison (log$_2$ fold change)'

    cm = sns.clustermap(
        df, row_colors=row_colors, cmap=cmap, norm=norm, clip_on=True,
        col_cluster=False, row_cluster=True, figsize=figsize,
        row_linkage=taxonomy_linkage(taxonomies), dendrogram_ratio=dendro_ratio,
    )
    for tick in cm.ax_row_colors.get_xticklabels():
        tick.set_fontsize(22)

    cbar = cm.ax_cbar
    fc_ticks = np.linspace(vmin, vmax, 5)
    cbar.set_yticks(fc_ticks)
    cbar.set_yticklabels([f'{t:+.1f}' for t in fc_ticks], fontsize=22)
    cbar.set_xlabel('log$_2$ fold change', fontsize=18, labelpad=10)

    cm.ax_col_dendrogram.set_visible(False)
    cm.figure.subplots_adjust(bottom=0.15, top=0.95)
    cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(), fontsize=24, rotation=0)
    cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(), fontsize=24, rotation=70, ha='right')
    cm.ax_heatmap.set_xlabel(xlabel, fontsize=32, labelpad=20)
    cm.ax_row_dendrogram.xaxis.set_visible(False)
    cm.ax_row_dendrogram.text(0.65, -0.02, 'Taxonomical tree', fontsize=24, ha='center',
                              transform=cm.ax_row_dendrogram.transAxes)

    # [phylum colors] [heatmap] [labels] [compressed dendrogram]
    colors_left = 0.18
    colors_w = 0.018
    heatmap_left = colors_left + colors_w + 0.005
    hm_pos = cm.ax_heatmap.get_position()
    dend_pos = cm.ax_row_dendrogram.get_position()
    cm.ax_heatmap.set_position([heatmap_left, hm_pos.y0, heatmap_right - heatmap_left, hm_pos.height])
    cm.ax_heatmap.yaxis.tick_right()
    cm.ax_heatmap.yaxis.set_label_position('right')
    cm.ax_heatmap.tick_params(axis='y', pad=4)
    cm.ax_row_dendrogram.set_position([heatmap_right + 0.10, dend_pos.y0, dendro_w, dend_pos.height])
    cm.ax_row_dendrogram.invert_xaxis()  # branches open leftward, toward heatmap

    col_dend_pos = cm.ax_col_dendrogram.get_position()
    cbar_w = 0.035
    cbar_x = colors_left - cbar_w - 0.04
    cm.ax_cbar.set_position([cbar_x, col_dend_pos.y0, cbar_w, col_dend_pos.height])

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

    # Move phylum color strip to the LEFT of the heatmap.
    hm_pos = cm.ax_heatmap.get_position()
    strip_w = 0.015
    cm.ax_row_colors.set_position([hm_pos.x0 - strip_w - 0.005, hm_pos.y0, strip_w, hm_pos.height])

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
                    handles.append(Patch(facecolor=color, label=f'      {cls}'))
            else:
                handles.append(Patch(facecolor=phylum_color[p], label=p))
    ncol = int(np.ceil(len(handles) / 2))
    # Sit the legend's bottom just above the top edge of the heatmap.
    hm_top = cm.ax_heatmap.get_position().y1
    cm.figure.legend(
        handles=handles, title='Phylum', title_fontsize=20, fontsize=16,
        loc='lower center', bbox_to_anchor=(0.5, hm_top + 0.005), ncol=ncol,
        frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6,
        columnspacing=1.2,
    )
    for spine in cm.ax_heatmap.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(1)

    if mode == 'log2fc_days' and phase_day_ranges:
        _draw_phase_delimiters(cm, df, phase_day_ranges)

    # Anchor the dendrogram's left edge to the rightmost edge of the longest label.
    fig = cm.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    max_right_disp = None
    for lab in cm.ax_heatmap.get_yticklabels():
        if not lab.get_text():
            continue
        bb = lab.get_window_extent(renderer=renderer)
        if max_right_disp is None or bb.x1 > max_right_disp:
            max_right_disp = bb.x1
    if max_right_disp is not None:
        max_right_fig = fig.transFigure.inverted().transform((max_right_disp, 0))[0]
        d_pos = cm.ax_row_dendrogram.get_position()
        cm.ax_row_dendrogram.set_position([max_right_fig + 0.008, d_pos.y0, dendro_w, d_pos.height])

    out_path = f"{title.lower().replace(' ', '_')}{suffix}.png"
    cm.figure.savefig(out_path, bbox_inches='tight', dpi=300)
    print(f'wrote {out_path}')


def main():
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

    # Aggregate to root organism ID and keep taxa with max relative abundance ≥ 1%.
    root = abundances.copy()
    root.columns = [c.split('.')[0] for c in root.columns]
    root = root.T.groupby(level=0).sum().T
    root = root.loc[:, root.max(axis=0) >= MIN_MAX_PCT]
    root_taxonomies = {}
    for ID, taxa in iterativeID_taxonomy.items():
        r = ID.split('.')[0]
        if r in root.columns and r not in root_taxonomies:
            root_taxonomies[r] = '|'.join(str(taxa[l]) for l in TAXONOMIC_LEVELS[:GENUS_DEPTH + 1])

    shared_kwargs = dict(
        iterativeID_color_map=iterativeID_color_map,
        genera_color_map=genera_color_map,
        proteo_class_color=proteo_class_color,
        proteo_base=proteo_base,
        iterativeID_levels=iterativeID_levels,
        phase_day_ranges=phase_day_ranges,
    )

    # ---- Phase-pair log2 fold change ----
    phase_means = {}
    for phase, span in phase_day_ranges.items():
        in_phase = [d for d in root.index if span['start'] <= int(d) <= span['end']]
        if in_phase:
            phase_means[phase] = root.loc[in_phase].mean(axis=0)
    phase_mean_df = pd.DataFrame(phase_means)
    phase_order = [p for p in phase_day_ranges if p in phase_mean_df.columns]
    fc_cols = {}
    for i, p1 in enumerate(phase_order):
        for p2 in phase_order[i + 1:]:
            fc_cols[f'{p1} / {p2}'] = np.log2(
                (phase_mean_df[p1] + PSEUDO) / (phase_mean_df[p2] + PSEUDO))
    df_phase = pd.DataFrame(fc_cols).astype(float).replace([np.inf, -np.inf], np.nan)
    tax_phase = pd.Series({idx: root_taxonomies.get(idx, f'Unknown|{idx}') for idx in df_phase.index})
    create_heatmap(df_phase, tax_phase, 'Taxa Above 1% Max Abundance',
                   mode='log2fc', suffix='_diff', **shared_kwargs)

    # ---- Per-day log2 fold change vs the first recorded day (innoculum) ----
    first_day = root.index[0]
    first_row = root.loc[first_day]
    df_innoc = pd.DataFrame({
        d: np.log2((root.loc[d] + PSEUDO) / (first_row + PSEUDO))
        for d in root.index
    }).astype(float).replace([np.inf, -np.inf], np.nan)
    tax_innoc = pd.Series({idx: root_taxonomies.get(idx, f'Unknown|{idx}') for idx in df_innoc.index})
    create_heatmap(df_innoc, tax_innoc, 'Taxa Above 1% Max Abundance',
                   mode='log2fc_days', suffix='_diff_innoculum', **shared_kwargs)


if __name__ == '__main__':
    main()
