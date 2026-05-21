"""Standalone driver that builds the co-occurrence network and renders the
_grid_layout variant by invoking the same render_network helper as the
co-occurrence cell in data_processing.py.

Run from the EmilyKin repo root with codiffusion_bioreactor's .venv:
    /home/andrew/repos/codiffusion_bioreactor/.venv/bin/python render_cooccurrence_grid.py
"""
import json
import math as _math
from collections import defaultdict
from itertools import combinations

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from numpy import ones as _np_ones
from pandas import DataFrame, read_csv
from scipy.spatial import ConvexHull as _ConvexHull
from scipy.stats import spearmanr


GAOs_PAOs = {
    'PAOs': ['Ca_Accumulibacter', 'Tetrasphaera', 'Dechloromonas', 'Microlunatus', 'Azonexus', 'Ca_Phosphoribacter'],
    'GAOs': ['Ca_Competibacter', 'Defluviicoccus', 'Propionivibrio', 'Ca_Contendobacter'],
    'Putative PAOs': ['Ca_Obscuribacter', 'Thauera', 'Zoogloea', 'Paracoccus'],
    'Putative GAOs': ['Micropruina', 'Amaricoccus', 'Ca_Glycocaulis', 'Thauera'],
    'Other PHA storing potential+ function': ['Pseudomonas', 'Bacillus', 'Acinetobacter', 'Rhodocyclaceae'],
}

iterativeID_taxonomy_local = json.load(open('iterativeID_taxonomy.json'))
level_3 = 'Phylum'
iterativeID_level = {_ID: content.get(level_3, 'Unknown') for _ID, content in iterativeID_taxonomy_local.items()}

_significantly_connected_organisms = [str(x) for x in np.load('FDR_passing_pairs.npy')]
_significantly_connected_organisms.append('Methanobacteriaceae.1')

iterativeID_color_map_local = json.load(open('iterativeID_color_map.json', 'r'))
order_color_map_local = json.load(open('Phylum_color_map.json', 'r'))

abund_thresh = 0.005
abundances_3 = read_csv('abundances.csv', header=0).set_index('sample')
abundances_3.drop([col for col in abundances_3.columns if col not in _significantly_connected_organisms], axis=1, inplace=True)


def _corr_pvalues_local(df):
    # vectorized: scipy.stats.spearmanr on whole matrix at once
    res = spearmanr(df.values, axis=0, nan_policy='omit')
    rho_mat = np.asarray(res.correlation)
    p_mat = np.asarray(res.pvalue)
    rho_mat = np.where(np.isnan(rho_mat), 0.0, rho_mat)
    p_mat = np.where(np.isnan(p_mat), 1.0, p_mat)
    cr = DataFrame(rho_mat, index=df.columns, columns=df.columns)
    pv = DataFrame(p_mat, index=df.columns, columns=df.columns)
    return pv, cr


pvals_matrix, corr_matrix_2 = _corr_pvalues_local(abundances_3)
_mean_rel_abund = abundances_3.mean(axis=0)

_presence = (abundances_3 > 0).astype(int)
_cooccurrence = defaultdict(int)
for _sample in _presence.itertuples(index=False):
    _present = [col for col, val in zip(_presence.columns, _sample) if val]
    for _pair in combinations(sorted(_present), 2):
        _cooccurrence[_pair] += 1

G = nx.Graph()
for (_a, _b), _count in _cooccurrence.items():
    _rho = corr_matrix_2.loc[_a, _b]
    _p_value = pvals_matrix.loc[_a, _b]
    if np.isnan(_rho):
        continue
    G.add_edge(_a, _b, weight=np.abs(_rho), rho=_rho, pvalue=_p_value, cooccurrence=_count)

print(f'Edges in G: {G.number_of_edges()}')
print(f'Nodes in G: {G.number_of_nodes()}')

_inverted_GAOs_PAOs = {_v: _k for _k, vs in GAOs_PAOs.items() for _v in vs}
_GAO_PAO_label_colors = {
    'GAOs': 'green',
    'Putative GAOs': 'mediumseagreen',
    'PAOs': 'blue',
    'Putative PAOs': 'cornflowerblue',
    'Other PHA storing potential+ function': 'red',
}


def _gao_pao_label_color(node_id):
    text = str(node_id)
    for _org, _cat in _inverted_GAOs_PAOs.items():
        if _org in text:
            return _GAO_PAO_label_colors.get(_cat, 'black')
    return 'black'


def header_patch(title):
    return mpatches.Patch(color='none', label=f'$\\bf{{{title}}}$')


def render_network(G_sub, title_slug, focus_nodes=None, grid_layout_mode=False, use_module_weighting=True, resolution=1.0, simplified=False):
    if G_sub.number_of_edges() == 0:
        print(f'[{title_slug}] no edges — skipping')
        return

    _pos_edges = [(u, v) for u, v, d in G_sub.edges(data=True) if d['rho'] > 0]
    _pos_subgraph = G_sub.edge_subgraph(_pos_edges).copy()
    try:
        _louvain_comms = nx.community.louvain_communities(
            _pos_subgraph,
            weight='weight', resolution=resolution, seed=42)
    except AttributeError:
        _louvain_comms = list(nx.algorithms.community.greedy_modularity_communities(G_sub, weight='weight'))
    _louvain_comms = list(_louvain_comms)
    try:
        _Q_for_label = nx.community.modularity(_pos_subgraph, _louvain_comms, weight='weight') if _louvain_comms else 0.0
    except AttributeError:
        _Q_for_label = nx.algorithms.community.modularity(_pos_subgraph, _louvain_comms, weight='weight') if _louvain_comms else 0.0

    _layout_comms = list(_louvain_comms)
    _node_module = {n: i for i, c in enumerate(_layout_comms) for n in c}
    _isolate_nodes = [n for n in G_sub.nodes() if n not in _node_module]
    for _iso in _isolate_nodes:
        _node_module[_iso] = len(_layout_comms)
        _layout_comms.append({_iso})
    if _isolate_nodes:
        print(f'[{title_slug}] isolates: {len(_isolate_nodes)} -> singleton modules')

    _multi_comms_only = [c for c in _layout_comms if len(c) >= 2]
    _n_singletons = sum(1 for c in _layout_comms if len(c) == 1)

    _module_export = {}
    _louvain_n = 0
    _singleton_n = 0
    for _comm in sorted(_layout_comms, key=lambda c: -len(c)):
        _members = sorted(_comm)
        if len(_comm) >= 2:
            _louvain_n += 1
            _module_export[f'module_{_louvain_n}'] = {'size': len(_comm), 'type': 'louvain', 'members': _members}
        else:
            _singleton_n += 1
            _module_export[f'singleton_{_singleton_n}'] = {'size': 1, 'type': 'isolate_singleton', 'members': _members}
    json.dump(_module_export, open(f'network_module_membership_{title_slug}.json', 'w'), indent=2)
    print(f'[{title_slug}] wrote network_module_membership_{title_slug}.json ({_louvain_n} louvain + {_singleton_n} singletons)')

    np.random.seed(42)
    _anchor_radius = 9.0
    _initial_pos = {}
    _multi_comms_sorted = [c for c in sorted(_layout_comms, key=len, reverse=True) if len(c) >= 2]
    for _i, _comm in enumerate(_multi_comms_sorted):
        _angle = 2 * np.pi * _i / max(len(_multi_comms_sorted), 1)
        _anchor = np.array([np.cos(_angle), np.sin(_angle)]) * _anchor_radius
        _jitter_scale = 0.30 + 0.02 * len(_comm)
        for _n in _comm:
            _initial_pos[_n] = _anchor + np.random.normal(0, _jitter_scale, 2)
    for _n in G_sub.nodes():
        if _n not in _initial_pos:
            _initial_pos[_n] = np.random.normal(0, 0.7, 2)

    if use_module_weighting:
        _module_layouts = {}
        _scale_factor = 35.0 if grid_layout_mode else 20.0
        for _mi, _comm in enumerate(_layout_comms):
            if len(_comm) < 2:
                continue
            _sub = G_sub.subgraph(_comm).copy()
            _module_scale = (len(_comm) if grid_layout_mode else _math.sqrt(len(_comm))) * _scale_factor
            if _sub.number_of_edges() == 0:
                _module_layouts[_mi] = {n: np.random.normal(0, _module_scale * 0.3, 2) for n in _comm}
            else:
                _module_layouts[_mi] = nx.spring_layout(
                    _sub, seed=42,
                    iterations=500 if grid_layout_mode else 300,
                    k=30.0 if grid_layout_mode else 5.0,
                    scale=_module_scale,
                    weight=None if grid_layout_mode else 'weight')
        if grid_layout_mode:
            _intra_amp = 5.4 * 3
            for _mi in _module_layouts:
                for _n in list(_module_layouts[_mi].keys()):
                    _module_layouts[_mi][_n] = _module_layouts[_mi][_n] * _intra_amp
            _max_module_radius = max(
                (max(np.linalg.norm(p) for p in _module_layouts[_mi].values())
                 for _mi in _module_layouts), default=1.0)
            _hull_buffer = 1.4
            _gap_factor = 1.0 / 1.05
            _n_multi_mods = len(_module_layouts)
            _grid_n = int(_math.ceil(_math.sqrt(max(_n_multi_mods, 1))))
            _cell_size = 2 * _max_module_radius * _hull_buffer * _gap_factor
            _v_stretch = 1.6
            _center_left_shift = -0.85 * _cell_size
            _right_v_extra = 1.4
            print(f'[{title_slug}] rotated grid: cell={_cell_size:.2f}, max module radius={_max_module_radius:.2f}, intra_amp={_intra_amp}')
            pos = {}
            _sorted_mis = sorted(_module_layouts.keys(), key=lambda mi: -len(_layout_comms[mi]))
            for _i, _mi in enumerate(_sorted_mis):
                _pre_row = _i // _grid_n
                _pre_col = _i % _grid_n
                _pre_col = {0: 1, 1: 0}.get(_pre_col, _pre_col)
                _pre_x = (_pre_col - (_grid_n - 1) / 2.0) * _cell_size
                _pre_y = (_pre_row - (_grid_n - 1) / 2.0) * _cell_size
                _cx = -_pre_y
                _cy = _pre_x * _v_stretch
                if abs(_cx) < 0.5 * _cell_size:
                    _cx += _center_left_shift
                elif _cx > 0.5 * _cell_size:
                    _cy *= _right_v_extra
                _module_center = np.array([_cx, _cy])
                for _n, _local in _module_layouts[_mi].items():
                    pos[_n] = _module_center + np.array([-_local[1], _local[0]])
            _singleton_mis = [_mi for _mi, _c in enumerate(_layout_comms) if len(_c) == 1]
            if _singleton_mis:
                _xext = _cell_size + _max_module_radius
                _yext = (_v_stretch * _right_v_extra) * _cell_size + _max_module_radius
                _singleton_radius = max(_xext, _yext) + _max_module_radius * 0.05
                for _i, _mi in enumerate(_singleton_mis):
                    _theta = 2 * np.pi * (_i + 0.5) / max(len(_singleton_mis), 1)
                    for _n in _layout_comms[_mi]:
                        pos[_n] = np.array([np.cos(_theta), np.sin(_theta)]) * _singleton_radius
        else:
            _max_module_radius = max(
                (max(np.linalg.norm(p) for p in _module_layouts[_mi].values())
                 for _mi in _module_layouts), default=1.0)
            _hull_buffer = 1.02
            _gap_factor = 1.01
            _n_multi_mods = len(_module_layouts)
            _min_chord = 2 * _max_module_radius * _hull_buffer * _gap_factor
            _sep_radius = _min_chord / (2 * np.sin(np.pi / _n_multi_mods)) if _n_multi_mods > 1 else 0
            print(f'[{title_slug}] ring radius={_sep_radius:.2f}')
            pos = {}
            _sorted_mis = sorted(_module_layouts.keys(), key=lambda mi: -len(_layout_comms[mi]))
            for _i, _mi in enumerate(_sorted_mis):
                _theta = 2 * np.pi * _i / max(_n_multi_mods, 1)
                _module_center = np.array([np.cos(_theta), np.sin(_theta)]) * _sep_radius
                for _n, _local in _module_layouts[_mi].items():
                    pos[_n] = _module_center + _local
            _singleton_mis = [_mi for _mi, _c in enumerate(_layout_comms) if len(_c) == 1]
            if _singleton_mis:
                _singleton_radius = _sep_radius + _max_module_radius * 0.3
                for _i, _mi in enumerate(_singleton_mis):
                    _theta = 2 * np.pi * (_i + 0.5) / max(len(_singleton_mis), 1)
                    for _n in _layout_comms[_mi]:
                        pos[_n] = np.array([np.cos(_theta), np.sin(_theta)]) * _singleton_radius
    else:
        pos = nx.spring_layout(G_sub, pos=_initial_pos, seed=42, iterations=500,
                               k=6.0, scale=7.5, weight='weight')

    norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = cm.RdBu
    width, height = (40, 30)
    fig, ax = plt.subplots(figsize=(width, height))
    scale = 5000 * (width / 10) * 2
    compressor = np.sqrt

    # --- Nodes sorted by mean rel. abundance ascending: most abundant drawn last (on top) ---
    _nodes_sorted = sorted(G_sub.nodes(), key=lambda n: _mean_rel_abund.get(n, 0))
    node_sizes = [scale * compressor(_mean_rel_abund.get(n, 0)) for n in _nodes_sorted]
    node_colors = [order_color_map_local.get(iterativeID_level.get(n, 'Unknown'), 'lightgray') for n in _nodes_sorted]
    edgecolors = [_gao_pao_label_color(n) if focus_nodes and n in focus_nodes else 'none' for n in _nodes_sorted]
    linewidths = [6 if focus_nodes and n in focus_nodes else 0 for n in _nodes_sorted]
    nx.draw_networkx_nodes(G_sub, pos, nodelist=_nodes_sorted, node_size=node_sizes, node_color=node_colors, edgecolors=edgecolors, linewidths=linewidths, alpha=0.9, ax=ax)

    if simplified:
        # --- Two aggregate edges per module-pair (one for rho>0, one for rho<0) ---
        # thickness ∝ |Σρ| within each bucket, color = mean ρ within that bucket
        from collections import defaultdict as _dd
        _module_centroid = {}
        for _mi, _comm in enumerate(_layout_comms):
            _pts = np.array([pos[_n] for _n in _comm if _n in pos])
            if len(_pts) == 0:
                continue
            _module_centroid[_mi] = _pts.mean(axis=0)
        _pair_pos_sum = _dd(float)
        _pair_pos_count = _dd(int)
        _pair_neg_sum = _dd(float)
        _pair_neg_count = _dd(int)
        for u, v, d in G_sub.edges(data=True):
            _mu = _node_module.get(u)
            _mv = _node_module.get(v)
            if _mu is None or _mv is None or _mu == _mv:
                continue
            _key = (min(_mu, _mv), max(_mu, _mv))
            _r = d['rho']
            if _r > 0:
                _pair_pos_sum[_key] += _r
                _pair_pos_count[_key] += 1
            elif _r < 0:
                _pair_neg_sum[_key] += _r
                _pair_neg_count[_key] += 1
        _all_counts = list(_pair_pos_count.values()) + list(_pair_neg_count.values())
        _max_count = max(_all_counts, default=1)
        _max_log = np.log10(_max_count + 1)
        _width_scale = 30.0 / _max_log if _max_log > 0 else 30.0
        _offset_dist = 0.012 * max(np.linalg.norm(np.array(list(_module_centroid.values())).max(0) - np.array(list(_module_centroid.values())).min(0)), 1.0)

        def _draw_split_edge(_key, _sum_rho, _count, _sign):
            if _key[0] not in _module_centroid or _key[1] not in _module_centroid:
                return
            _c_i = _module_centroid[_key[0]]
            _c_j = _module_centroid[_key[1]]
            _v = _c_j - _c_i
            _norm_v = np.linalg.norm(_v)
            if _norm_v < 1e-9:
                return
            _perp = np.array([-_v[1], _v[0]]) / _norm_v * _offset_dist * _sign
            _mean_rho = _sum_rho / _count
            ax.plot([_c_i[0] + _perp[0], _c_j[0] + _perp[0]],
                    [_c_i[1] + _perp[1], _c_j[1] + _perp[1]],
                    color=cmap(norm(_mean_rho)),
                    linewidth=_width_scale * np.log10(_count + 1),
                    alpha=0.85, solid_capstyle='round', zorder=2)

        # Draw smaller buckets first so larger ones overlay
        _all_buckets = (
            [(k, s, _pair_pos_count[k], +1) for k, s in _pair_pos_sum.items()]
            + [(k, s, _pair_neg_count[k], -1) for k, s in _pair_neg_sum.items()]
        )
        _all_buckets.sort(key=lambda x: x[2])
        for _key, _sum_rho, _count, _sign in _all_buckets:
            _draw_split_edge(_key, _sum_rho, _count, _sign)
        print(f'[{title_slug}] simplified: {len(_pair_pos_sum)} pos + {len(_pair_neg_sum)} neg inter-module edges (max count={_max_count})')
    else:
        # --- Edges sorted by |rho| ascending: stronger edges drawn last (on top) ---
        _edges_sorted = sorted(G_sub.edges(data=True), key=lambda e: abs(e[2]['rho']))
        _edgelist = [(u, v) for u, v, _ in _edges_sorted]
        rho_values = [d['rho'] for _, _, d in _edges_sorted]
        edge_widths = [3 * d['weight'] for _, _, d in _edges_sorted]
        edge_colors = [cmap(norm(r)) for r in rho_values]
        nx.draw_networkx_edges(G_sub, pos, edgelist=_edgelist, width=edge_widths, edge_color=edge_colors, alpha=0.85, ax=ax)

    # --- Labels: per-node zorder by abundance rank, so most abundant is on top ---
    for _rank, n in enumerate(_nodes_sorted):
        x, y = pos[n]
        abund = _mean_rel_abund.get(n, 0)
        font_size = 6 + 14 * compressor(abund) / compressor(_mean_rel_abund.max()) * (width / 10)
        font_size = max(5, min(font_size, 48))
        _lbl_color = _gao_pao_label_color(n)
        _txt = ax.text(x, y, str(n), fontsize=font_size, color=_lbl_color, fontweight='bold',
                       ha='center', va='center', zorder=5 + _rank)
        _txt.set_path_effects([path_effects.Stroke(linewidth=max(1.5, font_size / 6), foreground='white'), path_effects.Normal()])

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    _cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    _cbar.set_label('Spearman ρ', fontsize=10 * (width / 10))
    _cbar.ax.tick_params(labelsize=8 * (width / 10), length=8, width=2)

    legend_entries = [(0.002, '0.2%'), (0.02, '2%'), (0.20, '20%')]
    sizes_pt2 = [scale * compressor(_a) for _a, _ in legend_entries]
    diameters_pt = [2 * np.sqrt(_s / np.pi) for _s in sizes_pt2]
    breather_pt = 8
    offsets_pt = [0.0]
    for _i in range(1, len(legend_entries)):
        offsets_pt.append(offsets_pt[-1] + (diameters_pt[_i - 1] + diameters_pt[_i]) / 2 + breather_pt)
    fig_h_pts = fig.get_figheight() * 72
    top_y = 0.9
    circle_x = 0.85
    label_x = 0.88
    ax.text(circle_x, top_y + 0.025, 'Mean rel. abundance',
            transform=fig.transFigure, va='bottom', fontweight='bold',
            fontsize=6 * (width / 10), clip_on=False)
    for (_abund, _slabel), _s, _off in zip(legend_entries, sizes_pt2, offsets_pt):
        _y = top_y - _off / fig_h_pts
        ax.scatter([circle_x], [_y], s=_s, color='slategray', alpha=0.9,
                   transform=fig.transFigure, clip_on=False)
        ax.text(label_x, _y, _slabel,
                transform=fig.transFigure, va='center',
                fontsize=5 * (width / 10), clip_on=False)

    if simplified:
        # --- Edge-width legend anchored to bottom-left corner ---
        # Spacing chosen so bottom-most (thickest) line touches y=0 with the thickest line's
        # bottom edge just inside the figure.
        _count_legend_dy = 0.035
        _count_legend_top_y = 4 * _count_legend_dy + 0.01  # 0.15
        _line_x_start = 0.0
        _line_x_end = 0.06
        _count_label_x = 0.07
        ax.text(_line_x_start, _count_legend_top_y + 0.025,
                'Inter-module edges',
                transform=fig.transFigure, va='bottom', fontweight='bold',
                fontsize=6 * (width / 10), clip_on=False)
        ax.text(_line_x_start, _count_legend_top_y + 0.005,
                '(aggregated count)',
                transform=fig.transFigure, va='bottom',
                fontsize=5 * (width / 10), clip_on=False)
        _legend_counts = [c for c in [1, 10, 100, 1000, 10000] if c <= _max_count * 1.5]
        for _i, _count in enumerate(_legend_counts):
            _y = _count_legend_top_y - _i * _count_legend_dy
            _lw = _width_scale * np.log10(_count + 1)
            ax.plot([_line_x_start, _line_x_end], [_y, _y],
                    transform=fig.transFigure, color='gray',
                    linewidth=_lw, alpha=0.85,
                    solid_capstyle='round', clip_on=False)
            ax.text(_count_label_x, _y, f'{_count:,}',
                    transform=fig.transFigure, va='center',
                    fontsize=5 * (width / 10), clip_on=False)

    phyla_in_sub = sorted({iterativeID_level.get(n, 'Unknown') for n in G_sub.nodes()})
    _archaea_markers = ('archae', 'methano', 'halobac')
    _archaea_in_sub = [p for p in phyla_in_sub if any(m in p.lower() for m in _archaea_markers)]
    _bacteria_in_sub = [p for p in phyla_in_sub if p not in _archaea_in_sub and p in order_color_map_local]
    archaea_patches = [mpatches.Patch(color=order_color_map_local[_p], label=_p) for _p in _archaea_in_sub if _p in order_color_map_local]
    bacteria_patches = [mpatches.Patch(color=order_color_map_local[_p], label=_p) for _p in _bacteria_in_sub]
    legend_handles = []
    if archaea_patches:
        legend_handles.append(header_patch('Archaea'))
        legend_handles.extend(archaea_patches)
    if bacteria_patches:
        legend_handles.append(header_patch('Bacteria'))
        legend_handles.extend(bacteria_patches)
    ax.legend(handles=legend_handles, title='Taxonomic ' + level_3, title_fontsize=8 * (width / 10),
              loc='upper left', bbox_to_anchor=(0.02, 0.88), bbox_transform=fig.transFigure,
              fontsize=7 * (width / 10), frameon=True)
    ax.axis('off')
    plt.tight_layout()

    _palette_multi = plt.cm.tab10(np.linspace(0, 1, 10))
    _module_palette = []
    _multi_idx = 0
    for _c in _layout_comms:
        if len(_c) >= 2:
            _module_palette.append(_palette_multi[_multi_idx % 10])
            _multi_idx += 1
        else:
            _module_palette.append((0.65, 0.65, 0.65, 1.0))

    def _convex_polys_intersect(p1, p2):
        for poly, other in ((p1, p2), (p2, p1)):
            for _i in range(len(poly)):
                _edge = poly[(_i + 1) % len(poly)] - poly[_i]
                _normal = np.array([-_edge[1], _edge[0]])
                _proj_a = poly @ _normal
                _proj_b = other @ _normal
                if _proj_a.max() < _proj_b.min() or _proj_b.max() < _proj_a.min():
                    return False
        return True

    _hull_data = {}
    _initial_scale = 2.2
    for _mi, _comm in enumerate(_layout_comms):
        if len(_comm) < 2:
            continue
        _pts = np.array([pos[_n] for _n in _comm if _n in pos])
        if len(_pts) < 2:
            continue
        _centroid = _pts.mean(axis=0)
        if len(_pts) >= 3:
            _h = _ConvexHull(_pts)
            _verts = _pts[_h.vertices]
        else:
            _r = max(0.05, np.linalg.norm(_pts - _centroid, axis=1).max())
            _theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
            _verts = _centroid + _r * np.column_stack([np.cos(_theta), np.sin(_theta)])
        _hull_data[_mi] = [_centroid, _verts, _initial_scale]

    _min_scale = 1.4
    for _it in range(120):
        _overlap_found = False
        _ids = list(_hull_data.keys())
        for _i_a in range(len(_ids)):
            _mi_a = _ids[_i_a]
            _ca, _va, _sa = _hull_data[_mi_a]
            _poly_a = _ca + (_va - _ca) * _sa
            for _i_b in range(_i_a + 1, len(_ids)):
                _mi_b = _ids[_i_b]
                _cb, _vb, _sb = _hull_data[_mi_b]
                _poly_b = _cb + (_vb - _cb) * _sb
                if _convex_polys_intersect(_poly_a, _poly_b):
                    _overlap_found = True
                    if _hull_data[_mi_a][2] > _min_scale:
                        _hull_data[_mi_a][2] = max(_min_scale, _hull_data[_mi_a][2] * 0.95)
                    if _hull_data[_mi_b][2] > _min_scale:
                        _hull_data[_mi_b][2] = max(_min_scale, _hull_data[_mi_b][2] * 0.95)
        if not _overlap_found:
            print(f'[{title_slug}] hull non-overlap converged after {_it + 1} iterations')
            break
    else:
        print(f'[{title_slug}] hull non-overlap floor reached at scale {_min_scale}')

    for _mi, (_ca, _verts, _scale_v) in _hull_data.items():
        _poly_pts = _ca + (_verts - _ca) * _scale_v
        _poly = mpatches.Polygon(_poly_pts, closed=True,
                                 facecolor=_module_palette[_mi],
                                 edgecolor=_module_palette[_mi],
                                 linewidth=2, alpha=0.22, zorder=0)
        ax.add_patch(_poly)

    _module_summary = (
        r'$\bf{Modularity}$' + '\n'
        + f'Q = {_Q_for_label:.3f}' + '\n'
        + f'modules = {len(_multi_comms_only)} louvain + {_n_singletons} singletons '
        + f'(total = {len(_layout_comms)})' + '\n'
        + f'sizes (louvain): {sorted([len(c) for c in _multi_comms_only], reverse=True)}'
    )
    ax.text(0.02, 0.97, _module_summary,
            transform=fig.transFigure,
            fontsize=8 * (width / 10),
            ha='left', va='top',
            bbox=dict(boxstyle='round,pad=0.5',
                      facecolor='white', edgecolor='black',
                      alpha=0.85, linewidth=1.5))

    _layout_suffix = '_grid_layout' if grid_layout_mode else ('' if use_module_weighting else '_no_module_weighting')
    _gamma_suffix = f'_gamma{resolution:g}'
    _simplified_suffix = '_simplified' if simplified else ''
    out = f'cooccurrence_network_{title_slug}_min{abund_thresh*100:g}pct{_layout_suffix}{_gamma_suffix}{_simplified_suffix}.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'[{title_slug}] wrote {out}')


if __name__ == '__main__':
    render_network(G, 'p_value_FDR', grid_layout_mode=True, resolution=1.115, simplified=True)
