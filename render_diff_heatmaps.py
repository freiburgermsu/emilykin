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
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent
os.chdir(REPO)

TAXONOMIC_LEVELS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
DEFAULT_COLOR = 'lightgray'
MIN_MAX_PCT = 0.01
PSEUDO = 1e-6
GENUS_DEPTH = TAXONOMIC_LEVELS.index('Genus')

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


def _draw_phase_delimiters(clustermap, df, phase_day_ranges, time_based=False):
    """Vertical dashed lines at the lower bound of each phase (skipping the first),
    plus a phase-name label at the horizontal midpoint of each region.

    time_based=True maps day → day-valued x coordinate (used when the heatmap is
    drawn with pcolormesh on the day axis); otherwise day → column-index x.
    """
    days_int = []
    for c in df.columns:
        try:
            days_int.append(int(c))
        except (ValueError, TypeError):
            return
    n_cols = len(days_int)

    if time_based:
        def day_to_x(d):
            if d < days_int[0]:
                return float(days_int[0])
            if d > days_int[-1]:
                return float(days_int[-1])
            return float(d)
    else:
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
    # Place each dashed line on the cell BOUNDARY between two adjacent phases
    # (right edge of the last cell in the earlier phase), rather than through
    # the middle of any cell.
    phase_items = list(phase_day_ranges.items())
    for i in range(len(phase_items) - 1):
        earlier_span = phase_items[i][1]
        later_span = phase_items[i + 1][1]
        days_earlier = [d for d in days_int if earlier_span['start'] <= d <= earlier_span['end']]
        days_later = [d for d in days_int if later_span['start'] <= d <= later_span['end']]
        if not days_earlier or not days_later:
            continue
        last_day = max(days_earlier)
        first_day = min(days_later)
        if time_based:
            # pcolormesh edge between adjacent samples is at their midpoint.
            x_data = (last_day + first_day) / 2.0
        else:
            # Column-index mode: each cell occupies [col, col+1] in data coords;
            # boundary is at col + 1 (the right edge of the last cell).
            try:
                x_data = days_int.index(last_day) + 1.0
            except ValueError:
                continue
        x_disp = clustermap.ax_heatmap.transData.transform((x_data, 0))[0]
        x_fig = fig.transFigure.inverted().transform((x_disp, 0))[0]
        fig.add_artist(mlines.Line2D(
            [x_fig, x_fig], [y_bot, y_top], transform=fig.transFigure,
            color='black', linewidth=1.2, linestyle='--', alpha=0.6))

    # Phase-name (Roman numeral) labels disabled per request.
    # for phase_name, span in phase_day_ranges.items():
    #     days_in_phase = [d for d in days_int if span['start'] <= d <= span['end']]
    #     if not days_in_phase:
    #         continue
    #     mid_day = (min(days_in_phase) + max(days_in_phase)) / 2
    #     center_x_data = day_to_x(mid_day)
    #     if center_x_data is None:
    #         continue
    #     center_x_disp = clustermap.ax_heatmap.transData.transform((center_x_data, 0))[0]
    #     center_x_fig = fig.transFigure.inverted().transform((center_x_disp, 0))[0]
    #     fig.text(
    #         center_x_fig, y_top - 0.006, phase_name,
    #         ha='center', va='top', fontsize=36, fontweight='bold',
    #         color='black', transform=fig.transFigure,
    #     )


def create_heatmap(df, taxonomies, title, *, mode, suffix,
                   iterativeID_color_map, genera_color_map,
                   proteo_class_color, proteo_base, iterativeID_levels,
                   phase_day_ranges, significance=None):
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

    is_day_mode = mode in ('log2fc_days', 'log2fc_time')
    if is_day_mode:
        figsize = (50, 24)
        dendro_ratio = (0.06, 0.15)
        heatmap_right = 0.72
        dendro_w = 0.04
        xlabel = 'Days of operation'
    else:
        figsize = (max(14, 1.5 * df.shape[1] + 8), 24)
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
        tick.set_fontsize(26)

    cbar = cm.ax_cbar
    fc_ticks = np.linspace(vmin, vmax, 5)
    cbar.set_yticks(fc_ticks)
    cbar.set_yticklabels([f'{t:+.1f}' for t in fc_ticks], fontsize=26)
    cbar.set_xlabel('')  # clear default bottom label; the vertical label goes on the left
    cbar.set_ylabel('log$_2$ fold change', rotation=90, fontsize=22, labelpad=15)
    cbar.yaxis.set_label_position('left')
    cbar.yaxis.tick_right()  # tick labels on the right of the bar (toward the colors strip)

    cm.ax_col_dendrogram.set_visible(False)
    cm.figure.subplots_adjust(bottom=0.15, top=0.95)
    cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_yticklabels(), fontsize=29, rotation=0)
    cm.ax_heatmap.set_xticklabels(
        cm.ax_heatmap.get_xticklabels(), fontsize=29, rotation=80,
        ha='right', va='center', rotation_mode='anchor',
    )
    cm.ax_heatmap.set_xlabel(xlabel, fontsize=38, labelpad=20)
    cm.ax_row_dendrogram.xaxis.set_visible(False)
    cm.ax_row_dendrogram.text(0.65, -0.02, 'Taxonomical tree', fontsize=29, ha='center',
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

    # Initial cbar placement (height halved); final x is re-anchored to the
    # right of the dendrogram after labels render.
    cbar_w = 0.015
    cbar_h = hm_pos.height * 0.525  # 50% taller than the previous 0.35
    cbar_x = colors_left - cbar_w - 0.05
    cbar_y = hm_pos.y0 + (hm_pos.height - cbar_h) / 2
    cm.ax_cbar.set_position([cbar_x, cbar_y, cbar_w, cbar_h])

    if mode == 'log2fc_time':
        # Re-render the heatmap with pcolormesh so each column's width is
        # proportional to the time delta to its neighbors. Column width / time
        # delta is constant across the x-axis: x positions ARE the day values.
        ax = cm.ax_heatmap
        days = np.array([int(c) for c in df.columns], dtype=float)
        n_cols_t = len(days)
        if n_cols_t > 1:
            edges = np.empty(n_cols_t + 1)
            edges[1:-1] = (days[:-1] + days[1:]) / 2.0
            edges[0] = days[0] - (days[1] - days[0]) / 2.0
            edges[-1] = days[-1] + (days[-1] - days[-2]) / 2.0
        else:
            edges = np.array([days[0] - 0.5, days[0] + 0.5])
        reordered_idx = cm.dendrogram_row.reordered_ind
        reordered_data = df.values[reordered_idx, :]
        masked = np.ma.masked_invalid(reordered_data)
        n_rows_t = len(reordered_data)
        y_edges = np.arange(n_rows_t + 1, dtype=float)
        for _img in list(ax.images):
            _img.remove()
        for _coll in list(ax.collections):
            _coll.remove()
        ax.pcolormesh(edges, y_edges, masked, cmap=cmap, norm=norm, shading='flat')
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(n_rows_t, 0)
        # X-axis treated as a continuous graphical axis: ticks at every
        # multiple of 25 across the day domain, independent of which days
        # actually have samples.
        tick_step = 25
        x_lo, x_hi = float(edges[0]), float(edges[-1])
        first_tick = int(np.ceil(x_lo / tick_step) * tick_step)
        last_tick = int(np.floor(x_hi / tick_step) * tick_step)
        tick_days = list(range(first_tick, last_tick + 1, tick_step))
        ax.set_xticks(tick_days)
        ax.set_xticklabels(
            [str(d) for d in tick_days], fontsize=29, rotation=80,
            ha='right', va='center', rotation_mode='anchor',
        )
        ax.set_yticks(np.arange(n_rows_t) + 0.5)
        reordered_labels = [df.index[i] for i in reordered_idx]
        ax.set_yticklabels(reordered_labels, fontsize=29, rotation=0)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.tick_params(axis='y', pad=4)

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
    show_kingdom_headers = bool(archaea) and bool(bacteria)
    handles = []
    if archaea:
        if show_kingdom_headers:
            handles.append(Patch(color='none', label=r'$\bf{Archaea}$'))
        handles += [Patch(facecolor=phylum_color[p], label=p) for p in archaea]
    if bacteria:
        if show_kingdom_headers:
            handles.append(Patch(color='none', label=r'$\bf{Bacteria}$'))
        for p in bacteria:
            if p == 'Proteobacteria' and proteo_class_color:
                handles.append(Patch(color='none', label=p))
                for cls, color in proteo_class_color.items():
                    cls_short = cls.replace('proteobacteria', '').replace('Proteobacteria', '')
                    handles.append(Patch(facecolor=color, label=f'      {cls_short}'))
            else:
                handles.append(Patch(facecolor=phylum_color[p], label=p))
    # Phylum legend handles are built here; the legend itself is created later
    # so it can be anchored above the (final) right-side cbar.
    legend_handles = handles
    for spine in cm.ax_heatmap.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(1)

    if is_day_mode and phase_day_ranges:
        _draw_phase_delimiters(cm, df, phase_day_ranges, time_based=(mode == 'log2fc_time'))

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

    # Move the cbar to the far right of the figure: align its left edge with
    # the dendrogram's right edge (plus a small gap), keeping the halved height.
    dendro_pos = cm.ax_row_dendrogram.get_position()
    cbar_pos = cm.ax_cbar.get_position()
    cbar_target_left = dendro_pos.x1 + 0.030  # shifted further right
    cbar_final_y = hm_pos.y0  # anchor cbar at the bottom of the heatmap
    cm.ax_cbar.set_position([cbar_target_left, cbar_final_y, cbar_pos.width, cbar_pos.height])

    # Phylum legend (single column), top edge flush with the top of the figure
    # and horizontally centered on the cbar so they remain vertically in-line.
    # Text size matches the heatmap axis tick labels (29 pt).
    cbar_final_pos = cm.ax_cbar.get_position()
    legend_anchor_x = cbar_final_pos.x0 + cbar_final_pos.width / 2
    cm.figure.legend(
        handles=legend_handles, title='Phylum', title_fontsize=36, fontsize=29,
        loc='upper center', bbox_to_anchor=(legend_anchor_x, 0.84),
        ncol=1, frameon=True, borderaxespad=0.5,
        handlelength=1.5, handletextpad=0.6, columnspacing=1.2,
    )

    # Significance asterisks: drawn at each phase's horizontal midpoint, in data
    # coordinates of the time-based heatmap, for each (phase, taxon) cell whose
    # FDR-adjusted q-value passed.
    if significance and mode == 'log2fc_time' and phase_day_ranges:
        reordered_idx = cm.dendrogram_row.reordered_ind
        taxon_to_row = {df.index[reordered_idx[i]]: i + 0.5 for i in range(len(reordered_idx))}
        days_in_df = [int(c) for c in df.columns]
        phase_centers = {}
        for phase_name, span in phase_day_ranges.items():
            in_phase = [d for d in days_in_df if span['start'] <= d <= span['end']]
            if not in_phase:
                continue
            phase_centers[phase_name] = (min(in_phase) + max(in_phase)) / 2.0
        ax = cm.ax_heatmap
        for (phase_name, taxon), _q in significance.items():
            if phase_name not in phase_centers or taxon not in taxon_to_row:
                continue
            ax.text(phase_centers[phase_name], taxon_to_row[taxon], '*',
                    color='lime', fontsize=34, fontweight='bold',
                    ha='center', va='center')

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
    root_full = abundances.copy()
    root_full.columns = [c.split('.')[0] for c in root_full.columns]
    root_full = root_full.T.groupby(level=0).sum().T
    root = root_full.loc[:, root_full.max(axis=0) >= MIN_MAX_PCT]
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
    create_heatmap(df_innoc, tax_innoc, 'Taxa Above 1% Max Abundance',
                   mode='log2fc_time', suffix='_diff_innoculum_time', **shared_kwargs)

    # ---- CLR → one-sample Wilcoxon (phase vs innoculum) → BH-FDR ----
    # CLR on the full root-aggregated composition (not just filtered taxa).
    log_full = np.log(root_full + PSEUDO)
    clr_full = log_full.sub(log_full.mean(axis=1), axis=0)
    clr = clr_full[root.columns]  # subset to filtered taxa, rows still days
    innoculum_day = root.index[0]
    innoc_phase = next(
        (p for p, s in phase_day_ranges.items()
         if s['start'] <= int(innoculum_day) <= s['end']),
        None,
    )
    innoc_clr = clr.loc[innoculum_day]
    pvals_map = {}
    for phase, span in phase_day_ranges.items():
        if phase == innoc_phase:
            continue
        in_phase = [d for d in root.index if span['start'] <= int(d) <= span['end']]
        if len(in_phase) < 1:
            continue
        phase_clr = clr.loc[in_phase]
        for taxon in clr.columns:
            diffs = (phase_clr[taxon] - innoc_clr[taxon]).values
            nz = diffs[diffs != 0]
            if nz.size < 1:
                continue
            try:
                _, p = wilcoxon(nz, alternative='two-sided', zero_method='wilcox')
            except ValueError:
                continue
            pvals_map[(phase, taxon)] = p
    significance = {}
    if pvals_map:
        keys = list(pvals_map.keys())
        pvals = np.array([pvals_map[k] for k in keys])
        reject, qvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        significance = {keys[i]: qvals[i] for i in range(len(keys)) if reject[i]}
    print(f'FDR-significant (phase, taxon) cells: {len(significance)} / {len(pvals_map)}')
    create_heatmap(df_innoc, tax_innoc, 'Taxa Above 1% Max Abundance',
                   mode='log2fc_time', suffix='_diff_innoculum_time_FDR',
                   significance=significance, **shared_kwargs)


if __name__ == '__main__':
    main()
