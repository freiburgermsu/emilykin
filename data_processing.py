import marimo

__generated_with = "0.23.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Taxonomy and Abundance data
    """)
    return


@app.cell
def _(display):
    from json import load, dump
    from pandas import read_csv
    _total_df = read_csv('table_rel_export.csv').set_index('seq')
    _total_df = _total_df.fillna('')
    display(_total_df)
    relative_abundance, taxonomy, _sample_days = ({}, {}, {})
    for _seq, _row in _total_df.iterrows():
        relative_abundance.setdefault(_seq, {})
        taxonomy.setdefault(_seq, {'Kingdom': _row['Kingdom'], 'Phylum': _row['Phylum'], 'Class': _row['Class'], 'Order': _row['Order'], 'Family': _row['Family'], 'Genus': _row['Genus'], 'Species': _row['Species']})
        relative_abundance[_seq].setdefault(_row['sample'], 0)
        relative_abundance[_seq][_row['sample']] = relative_abundance[_seq][_row['sample']] + _row['rel_ab']
        _sample_days.setdefault(_row['sample'], _row['day'])
    dump(relative_abundance, open('relative_abundance.json', 'w'))
    dump(taxonomy, open('taxonomy.json', 'w'))
    dump(_sample_days, open('sample_days.json', 'w'))
    return dump, load, read_csv, taxonomy


@app.cell
def _(display, dump, taxonomy):
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    iterativeIDs, iterativeID_levels, _iterativeID_phylums, taxonIDs, _iterativeID_taxonomy = ({}, {}, {}, {}, {})
    for _seq, _taxa in taxonomy.items():
        best_level = 'Kingdom'
        for level, taxon in _taxa.items():
            if taxon == '' or level == 'Species':
                break
            best_level = level
        OGtaxon = _taxa[best_level]
        taxon = OGtaxon + '.1'
        while taxon in iterativeIDs:
            taxon, _count = ('.'.join(taxon.split('.')[:-1]), int(taxon.split('.')[-1]))
            taxon = f'{taxon}.{_count + 1}'
        taxonIDs.setdefault(OGtaxon, []).append(_seq)
        iterativeIDs[taxon] = _seq
        iterativeID_levels[taxon] = best_level
        _iterativeID_phylums[taxon] = _taxa['Phylum']
        _iterativeID_taxonomy[taxon] = taxonomy.get(_seq, 'Unknown')
    dump(iterativeIDs, open('iterativeIDs.json', 'w'))
    dump(iterativeID_levels, open('iterativeID_levels.json', 'w'))
    dump(_iterativeID_phylums, open('iterativeID_phylums.json', 'w'))
    dump(taxonIDs, open('taxonIDs.json', 'w'))
    dump(_iterativeID_taxonomy, open('iterativeID_taxonomy.json', 'w'))
    print(len(iterativeIDs), 'unique orgs')
    print(len(taxonIDs), 'unique taxa')
    display(list(iterativeIDs.items())[-10:])
    display(list(iterativeID_levels.items())[-10:])
    display(list(_iterativeID_phylums.items())[-10:])
    return iterativeID_levels, iterativeIDs, level


@app.cell
def _(iterativeID_levels):
    sorted({x.split(".")[0] for x in iterativeID_levels.keys()})
    return


@app.cell
def _(display, iterativeIDs, load):
    from pandas import DataFrame
    _ab = load(open('relative_abundance.json', 'r'))
    ab_df = DataFrame(_ab) / 100
    _seqIDs = {_seq: _ID for _ID, _seq in iterativeIDs.items()}
    ab_df.columns = [_seqIDs[col] for col in ab_df.columns]
    print(ab_df.shape, min(ab_df.min()), max(ab_df.max()))
    not_normalized = [x for x, val in ab_df.T.sum().items() if round(val, 1) != 1]
    display(len(not_normalized), not_normalized)
    _zero_level = 1e-05
    ab_df = ab_df.loc[:, (ab_df.fillna(0) > _zero_level).sum() >= 6]
    display(ab_df.shape)
    display(ab_df.head())
    ab_df.index.name = 'sample'
    ab_df.to_csv('abundances.csv')
    return (DataFrame,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Abundance heatmaps
    """)
    return


@app.cell
def _(array, defaultdict, np):
    def taxonomy_linkage(taxonomy_series):
        """
        Build a linkage matrix that EXACTLY follows taxonomy hierarchy.
        Auto-detects depth from the taxonomy strings.
        """
        parsed = _taxonomy_series.fillna('unknown').astype(str).str.split('[;|,]', regex=True).apply(lambda lst: [x.strip() for x in lst or [] if x.strip()])
        ranks_n = max((len(_p) for _p in parsed))
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
                    next_cluster_id = next_cluster_id + 1
                    linkage_rows.append([cluster_id, idx_to_pos[idx], merge_distance, total_count + 1])
                    cluster_id = new_id
                    total_count = total_count + 1
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
                next_cluster_id = next_cluster_id + 1
                linkage_rows.append([cluster_id, sub_id, merge_distance, total_count + sub_count])
                cluster_id = new_id
                total_count = total_count + sub_count
            return (cluster_id, total_count)
        build_subtree(list(parsed.index), depth=0)
        Z = array(linkage_rows, dtype=float)
        if len(Z) > 0:
            for _i in range(1, len(Z)):
                if Z[_i, 2] < Z[_i - 1, 2]:
                    Z[_i, 2] = Z[_i - 1, 2]
        return Z
    GAOs_PAOs = {'PAOs': ['Ca_Accumulibacter', 'Tetrasphaera', 'Dechloromonas', 'Microlunatus', 'Azonexus', 'Ca_Phosphoribacter'], 'GAOs': ['Ca_Competibacter', 'Defluviicoccus', 'Propionivibrio', 'Ca_Contendobacter'], 'Putative PAOs': ['Ca_Obscuribacter', 'Thauera', 'Zoogloea', 'Paracoccus'], 'Putative GAOs': ['Micropruina', 'Amaricoccus', 'Ca_Glycocaulis', 'Thauera'], 'Other PHA storing potential+ function': ['Pseudomonas', 'Bacillus', 'Acinetobacter', 'Rhodocyclaceae']}
    return GAOs_PAOs, taxonomy_linkage


@app.cell
def _(dump, iterativeIDs, level, load, read_csv, taxonomy):
    from pathlib import Path
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches
    from collections import defaultdict
    from itertools import combinations
    import numpy as np
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import multipletests

    def rgba_to_hex(rgba):
        return '#{:02x}{:02x}{:02x}'.format(*(int(c * 255) for c in rgba[:3]))

    def build_phylum_colors(out_path='phylum_colors.json'):
        phyla = sorted({t.get('Phylum') for t in taxonomy.values() if t.get('Phylum')})
        archaea_markers = ['archaeo', 'euryarchaeota', 'crenarchaeota', 'thaumarchaeota', 'halobacterota', 'methanobacteriota', 'micrarchaeota', 'nanoarchaeota']
        archaea = [_p for _p in phyla if any((m in _p.lower() for m in archaea_markers))]
        bacteria = [_p for _p in phyla if _p not in archaea]
        colors = {}
        for _i, _p in enumerate(archaea):
            colors[_p] = rgba_to_hex(plt.cm.turbo(_i / max(len(archaea), 1) * 0.15))
        for _i, _p in enumerate(bacteria):
            colors[_p] = rgba_to_hex(plt.cm.turbo(0.2 + _i / max(len(bacteria), 1) * 0.8))
        dump(colors, open(out_path, 'w'), indent=2)
        return colors
    build_phylum_colors()
    _iterativeID_taxonomy = load(open('iterativeID_taxonomy.json', 'r'))
    _iterativeID_phylums = load(open('iterativeID_phylums.json', 'r'))
    none_keys = [k for k, v in _iterativeID_phylums.items() if v is None]
    print('None taxonomy entries:', none_keys)
    archaea_phyla = sorted({v for k, v in _iterativeID_phylums.items() if v is not None and ('archaeo' in v.lower() or any((_a in v.lower() for _a in ['candidatus thermoplasmatota', 'halobacterota', 'methanobacteriota'])))})
    bacteria_phyla = sorted({v for k, v in _iterativeID_phylums.items() if v is not None and v not in archaea_phyla})
    n_archaea = len(archaea_phyla)
    archaea_colors = [plt.cm.turbo(_i / max(n_archaea, 1) * 0.15) for _i in range(n_archaea)]
    n_bacteria = len(bacteria_phyla)
    bacteria_colors = [plt.cm.turbo(0.2 + _i / max(n_bacteria, 1) * 0.8) for _i in range(n_bacteria)]
    taxa_color_map = {}
    for _phylum, _color in zip(archaea_phyla, archaea_colors):
        taxa_color_map[_phylum] = _color
    for _phylum, _color in zip(bacteria_phyla, bacteria_colors):
        taxa_color_map[_phylum] = _color
    dump(taxa_color_map, open(f'{level}_color_map.json', 'w'))
    _iterativeID_color_map = {_ID: taxa_color_map[_phylum] for _ID, _phylum in _iterativeID_phylums.items() if _phylum}
    dump(_iterativeID_color_map, open(f'iterativeID_color_map.json', 'w'))
    phylum_colors = load(open('phylum_colors.json', 'r'))
    DEFAULT = 'lightgray'
    df = read_csv('abundances.csv').set_index('sample')
    df = df.loc[:, (df.fillna(0) > 0).sum() >= 6]
    _seqIDs = {k: v for v, k in iterativeIDs.items()}
    df.columns = [_seqIDs.get(_ID, _ID) for _ID in df.columns]
    _mean_rel_abund = df.div(df.sum(axis=1), axis=0).mean(axis=0)
    _zero_level = 1e-05
    _presence = (df > _zero_level).astype(int)
    _cooccurrence = defaultdict(int)
    for _sample in _presence.itertuples(index=False):
        _present = [col for col, val in zip(_presence.columns, _sample) if val]
        for _pair in combinations(sorted(_present), 2):
            _cooccurrence[_pair] = _cooccurrence[_pair] + 1
    pair_data = []
    for (_a, _b), _count in _cooccurrence.items():
        _rho, _p_value = spearmanr(df[_a], df[_b])
        if np.isnan(_rho):
            continue
        pair_data.append((_a, _b, _rho, _p_value, _count))
    _pvals = np.array([t[3] for t in pair_data])
    reject, qvals, _, _ = multipletests(_pvals, alpha=0.05, method='fdr_bh')
    passing = [pair_data[_i] for _i in range(len(pair_data)) if reject[_i]]
    mems = set()
    for _a, _b, _rho, _p_value, _count in passing:
        mems.add(_a)
        mems.add(_b)
    np.save('FDR_passing_pairs', list(mems))
    member_colors = {_iterativeID_phylums.get(n, n): _iterativeID_color_map.get(n, DEFAULT) for n in mems}
    archaea_phyla = {k for k in archaea_phyla if k in member_colors.keys()}
    bacteria_phyla = {k for k in bacteria_phyla if k in member_colors.keys()}
    archaea_patches = [mpatches.Patch(color=member_colors[_p], label=_p) for _p in archaea_phyla]
    bacteria_patches = [mpatches.Patch(color=member_colors[_p], label=_p) for _p in bacteria_phyla]
    return (
        archaea_phyla,
        bacteria_phyla,
        combinations,
        defaultdict,
        mcolors,
        mpatches,
        multipletests,
        np,
        pair_data,
        plt,
        reject,
        spearmanr,
    )


@app.cell
def _(
    DataFrame,
    GAOs_PAOs,
    display,
    genera_color_map,
    mcolors,
    np,
    read_csv,
    taxonomy_linkage,
):
    from pandas import Series, concat
    from matplotlib.colors import LogNorm, Normalize, BoundaryNorm, TwoSlopeNorm
    from matplotlib.patches import Patch
    import seaborn as sns
    from matplotlib import pyplot
    import json
    from numpy import inf, nan, empty, logspace, array, log10, nanmean, where, isnan
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage
    from collections import Counter
    import colorsys
    import sigfig
    _sample_days = json.load(open('sample_days.json'))
    _iterativeID_color_map = json.load(open('iterativeID_color_map.json'))
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    _abundances = read_csv('abundances.csv').set_index('sample')
    _abundances.index = [_sample_days[col] for col in _abundances.index]
    _abundances = _abundances.loc[sorted(_abundances.index, key=int)]
    intervals = [0] + list(logspace(-2, -0.3, 5))
    print(intervals)

    def create_heatmap(df, taxonomies, title, inlayed_data=None, linkage=None):
        new_cmap = LinearSegmentedColormap.from_list('NewMap', [(0.0, 'aliceblue'), (0.25, 'lightblue'), (1.0, 'navy')])
        new_cmap.set_bad('aliceblue')
        vmin = df.min().min()
        vmax = df.max().max()
        vcenter = log10(0.1)
        vcenter = min(max(vcenter, vmin + 1e-06), vmax - 1e-06)
        norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        DEFAULT_COLOR = 'lightgray'

        def _lighten(c, factor):
            h, l, s = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
            return colorsys.hls_to_rgb(h, max(0.0, min(1.0, l * factor)), s)
        proteo_base = _iterativeID_color_map.get('Proteobacteria', 'tab:purple')
        proteo_classes = sorted({str(_taxonomies.get(_i, '')).split('|')[2] for _i in df.index if len(str(_taxonomies.get(_i, '')).split('|')) >= 3 and str(_taxonomies[_i]).split('|')[1] == 'Proteobacteria'})
        n = max(len(proteo_classes), 1)
        proteo_class_color = {cls: _lighten(proteo_base, 0.6 + 0.8 * _i / max(n - 1, 1)) for _i, cls in enumerate(proteo_classes)}

        def lookup(idx):
            parts = str(_taxonomies.get(idx, '')).split('|')
            if len(parts) >= 3 and parts[1] == 'Proteobacteria':
                return proteo_class_color.get(parts[2], proteo_base)
            if idx in _iterativeID_color_map:
                return _iterativeID_color_map[idx]
            elif idx in genera_color_map:
                return genera_color_map[idx]
            return DEFAULT_COLOR
        row_colors = Series({idx: lookup(idx) for idx in df.index}, name='Phylum')
        _clusterMap = sns.clustermap(df, row_colors=row_colors, cmap=new_cmap, norm=norm, clip_on=True, col_cluster=False, row_cluster=True, figsize=(50, 20), row_linkage=taxonomy_linkage(_taxonomies), dendrogram_ratio=(0.2, 0.15))
        for tick in _clusterMap.ax_row_colors.get_xticklabels():
            tick.set_fontsize(22)
        _cbar = _clusterMap.ax_cbar
        log_ticks = np.array([vmin, log10(0.01), vcenter, log10(0.2), vmax])
        _cbar.set_yticks(log_ticks)
        original_ticks = 10 ** log_ticks
        original_labels = [f'{x * 100:.1e}'.replace('e-0', 'E-') if x * 100 < 0.1 else round(x * 100) for x in original_ticks]
        original_labels[0] = '0'
        _cbar.set_yticklabels(original_labels, fontsize=22)
        _cbar.set_yticks(log_ticks)
        _cbar.set_xlabel('Rel. Abundance %', fontsize=18, labelpad=10)
        if inlayed_data:
            top_ax = _clusterMap.ax_col_dendrogram
            top_ax.clear()
            top_ax.plot(inlayed_data.keys(), inlayed_data.values(), color='black', marker='o', linewidth=3)
            top_ax.set_ylabel('Shannon Diversity', fontsize=20)
            top_ax.yaxis.tick_right()
            top_ax.yaxis.set_label_position('right')
            top_ax.grid(True, axis='y', linestyle='--', alpha=0.7)
            top_ax.set_xticks([])
            top_ax.tick_params(axis='y', labelsize=16)
        _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
        _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=24, rotation=0)
        _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=24, rotation=70, ha='right')
        _clusterMap.ax_heatmap.set_xlabel('Days of operation', fontsize=32, labelpad=20)
        _clusterMap.ax_row_dendrogram.xaxis.set_visible(False)
        _clusterMap.ax_row_dendrogram.text(0.65, -0.02, 'Taxonomical tree', fontsize=24, ha='center', transform=_clusterMap.ax_row_dendrogram.transAxes)
        _clusterMap.ax_heatmap.yaxis.tick_left()
        _clusterMap.ax_heatmap.yaxis.set_label_position('left')
        label_right = 0.22
        heatmap_right = 0.68
        dendro_left = 0.72
        dendro_right = 0.85
        hm_pos = _clusterMap.ax_heatmap.get_position()
        dend_pos = _clusterMap.ax_row_dendrogram.get_position()
        _clusterMap.ax_cbar.set_position([0.09, 0.86, 0.06, 0.14])
        _clusterMap.ax_heatmap.set_position([label_right, hm_pos.y0, heatmap_right - label_right, hm_pos.height])
        _clusterMap.ax_row_dendrogram.set_position([dendro_left, dend_pos.y0, dendro_right - dendro_left, dend_pos.height])
        _clusterMap.ax_row_dendrogram.invert_xaxis()
        gap = 0.03
        col_dend_pos = _clusterMap.ax_col_dendrogram.get_position()
        _clusterMap.ax_col_dendrogram.set_position([label_right, col_dend_pos.y0 + gap, heatmap_right - label_right, col_dend_pos.height])
        iterativeID_levels = json.load(open('iterativeID_levels.json', 'r'))
        _ID_levels = {}
        for k, v in iterativeID_levels.items():
            _ID_levels[k] = v
            _ID_levels.setdefault(k.split('.')[0], v)
        inverted_GAOs_PAOs = {v: k for k, vs in GAOs_PAOs.items() for v in vs}
        for _label in _clusterMap.ax_heatmap.get_yticklabels():
            _text = _label.get_text()
            for _org, val in inverted_GAOs_PAOs.items():
                if _org not in _text:
                    continue
                _label.set_fontweight('bold')
                if 'GAOs' in val:
                    _label.set_color('green')
                    if 'Putative' in val:
                        _label.set_color('lightgreen')
                elif 'PAOs' in val:
                    _label.set_color('blue')
                    if 'Putative' in val:
                        _label.set_color('lightblue')
                elif 'PHA' in val:
                    _label.set_color('red')
            if _ID_levels.get(_text) == 'Genus':
                _label.set_fontstyle('italic')
        hm_pos = _clusterMap.ax_heatmap.get_position()
        fig_w = _clusterMap.figure.get_figwidth()
        fig_h = _clusterMap.figure.get_figheight()
        strip_w = 0.015
        strip_h = 0.015
        _clusterMap.ax_row_colors.set_position([hm_pos.x1 + 0.005, hm_pos.y0, strip_w, hm_pos.height])
        y_pad_pts = strip_w * fig_w * 72 + 12
        _clusterMap.ax_heatmap.tick_params(axis='y', pad=y_pad_pts)
        if isinstance(row_colors, Series):
            rc = row_colors
        else:
            rc = Series(list(row_colors), index=df.index)
        phylum_color = {}
        for idx in df.index:
            parts = str(_taxonomies.get(idx, '')).split('|')
            if len(parts) < 2:
                continue
            _phylum = parts[1]
            if _phylum in ('None', '', 'Unknown'):
                continue
            if _phylum == 'Proteobacteria':
                phylum_color['Proteobacteria'] = proteo_base
                continue
            _color = rc.get(idx)
            if _color is None:
                continue
            phylum_color.setdefault(_phylum, _color)
        print('phylum_color', phylum_color)
        archaea_markers = ('archaeo', 'halobacterota', 'methanobacteriota')

        def is_archaea(p):
            return any((m in _p.lower() for m in archaea_markers))
        archaea = sorted((_p for _p in phylum_color if is_archaea(_p)))
        bacteria = sorted((_p for _p in phylum_color if not is_archaea(_p)))
        handles = []
        if archaea:
            handles.append(Patch(color='none', label='$\\bf{Archaea}$'))
            handles = handles + [Patch(facecolor=phylum_color[_p], label=_p) for _p in archaea]
        if bacteria:
            handles.append(Patch(color='none', label='$\\bf{Bacteria}$'))
            for _p in bacteria:
                if _p == 'Proteobacteria' and proteo_class_color:
                    handles.append(Patch(color='none', label=_p))
                    for cls, _color in proteo_class_color.items():
                        handles.append(Patch(facecolor=_color, label=f'      {cls}'))
                else:
                    handles.append(Patch(facecolor=phylum_color[_p], label=_p))
        _clusterMap.figure.legend(handles=handles, title='Phylum', title_fontsize=20, fontsize=18, loc='upper right', bbox_to_anchor=(0.8, 1.0), frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6)
        for spine in _clusterMap.ax_heatmap.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor('black')
            spine.set_linewidth(1)
        _clusterMap.figure.savefig(f"{title.lower().replace(' ', '_')}.png", bbox_inches='tight', dpi=300)
    from numpy import log
    iterativeIDs_1 = json.load(open('iterativeIDs.json', 'r'))
    _zero_level = 1e-05

    def shannon_index(abundances):
        return -sum([abundance * log(abundance) for abundance in _abundances if abundance > _zero_level])
    shannon_indices = {_sample: shannon_index(list(abs.to_numpy())) for _sample, abs in _abundances.iterrows()}
    taxonomy_1 = json.load(open('iterativeID_taxonomy.json', 'r'))
    dic, _taxonomies = ({}, {})
    level_1 = 'Genus'
    for _sample, abs in _abundances.iterrows():
        for _org, _ab in abs.items():
            if _ab <= _zero_level:
                continue
            _taxonomies.setdefault(_org, '|'.join([str(taxonomy_1[_org][l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index('Genus')]))
            dic.setdefault(_sample, {})
            dic[_sample].setdefault(_org, 0)
            dic[_sample][_org] = dic[_sample][_org] + _ab
    _nonzero_per_day = {_day: dict(sorted({k: v for k, v in _org_dict.items() if v > 0}.items(), key=lambda item: item[1], reverse=True)) for _day, _org_dict in dic.items()}
    json.dump(_nonzero_per_day, open('nonzero_per_day.json', 'w'), indent=2)
    _top_per_day = {}
    _all_orgs = {}
    _topNum = 10
    top10_all_days = []
    for _day, _org_dict in _nonzero_per_day.items():
        top10_all_days.extend(list(_org_dict.keys())[:_topNum])
    top10_all_days = list(set(top10_all_days))
    for _day, _org_dict in _nonzero_per_day.items():
        _top_per_day.setdefault(_day, {})
        _top_per_day[_day] = {_org: log10(v) for _org, v in _org_dict.items() if _org in top10_all_days}
    _taxonomies = {_org: _taxa for _org, _taxa in _taxonomies.items() if _org in top10_all_days}
    df_1 = DataFrame(_top_per_day)
    df_1 = df_1.astype(float).replace([inf, -inf], nan)
    _taxonomy_series = Series({idx: _taxonomies.get(idx, f'Unknown|{idx}') for idx in df_1.index})
    display(df_1)
    _arr = df_1.values.astype(float)
    row_means = nanmean(_arr, axis=1, keepdims=True)
    arr_filled = where(isnan(_arr), row_means, _arr)
    row_linkage = linkage(arr_filled, method='average', metric='euclidean')
    create_heatmap(df_1, _taxonomy_series, f'Top {_topNum} ASVs (% abundance)', shannon_indices, row_linkage)
    return Patch, Series, abs, array, inf, json, nan, sns


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Correlate member abundances with methane
    """)
    return


@app.cell
def _(merged_normalized_abundances, normalized_methane):
    merged_correlations = {}
    for index, _row in merged_normalized_abundances.iterrows():
        corr = _row.corr(normalized_methane)
        merged_correlations[index] = corr
    merged_correlations = dict(sorted(merged_correlations.items(), key=lambda item: item[1], reverse=True))  # print(f"{index} correlates {corr} with methane")
    merged_correlations
    return


@app.cell
def _(DataFrame, Series, dump, load, spearmanr):
    from pandas import isna
    iterativeIDs_2 = load(open('model_inputs/iterativeIDs.json', 'r'))
    _abundances = load(open('model_inputs/abundances.json', 'r'))
    asvSet_abundances = load(open(f'modeling_files/ASVset_abundances.json', 'r'))
    omitted_columns = {'all_relevant_samples': ['10AB', 'A34', 'B34', 'C3', 'D34', 'E34', 'F34', 'G12', 'G3', 'H34']}

    def _correlations(ser1, ser2):
        aligned1, aligned2 = ser1.align(ser2, join='inner')
        if len(set(aligned1)) == 1 or len(set(aligned2)) == 1:
            global constant_vals
            constant_vals = constant_vals + 1
            return (float('nan'), float('nan'))
        return spearmanr(aligned1, aligned2)
    abundances_df = DataFrame(_abundances).drop(ommitted, axis=1)
    abundances_df = abundances_df.loc[:, (abundances_df.fillna(0) > 0).sum() >= 5]
    ASVset_abundances_df = DataFrame(asvSet_abundances).drop(ommitted, axis=1).fillna(0)
    ASVset_abundances_df = ASVset_abundances_df.loc[:, ASVset_abundances_df.notna().sum() >= 5]
    print(abundances_df.shape)
    for dataName, term in {'H2_BT': 'H2 breakthrough (mol/min)', 'H2_feed': 'H2 delivery rate (mol/min)'}.items():
        for _name, ommitted in omitted_columns.items():
            constant_vals = 0
            summary_samples = load(open('model_inputs/measurements/summary_samples.json', 'r'))
            sample_values = {_sample: content[term] for _sample, content in summary_samples.items() if _sample not in ommitted}
            print(sample_values)
            ASV_correlations = {}
            constant_vals = 0
            for ASV, sampleAbun in abundances_df.iterrows():
                _ID = iterativeIDs_2.get(ASV, ASV)
                correlation, _p = _correlations(sampleAbun, Series(sample_values))
                if isna(correlation):
                    continue
                ASV_correlations[_ID] = {'correlation': correlation, 'p_value': _p}
            ASV_correlations = dict(sorted(ASV_correlations.items(), key=lambda item: item[1]['correlation'], reverse=True))
            print(f'Constant values: {constant_vals}')
            dump(ASV_correlations, open(f'modeling_files/correlations/ASV_correlations_{_name}_{dataName}.json', 'w'))
            IterativeID_correlations = {}
            for ASV, sampleAbun in ASVset_abundances_df.iterrows():
                _ID = iterativeIDs_2.get(ASV, ASV)
                correlation, _p = _correlations(sampleAbun, Series(sample_values))
                if isna(correlation):
                    continue
                IterativeID_correlations[_ID] = {'correlation': correlation, 'p_value': _p}
            IterativeID_correlations = dict(sorted(IterativeID_correlations.items(), key=lambda item: item[1]['correlation'], reverse=True))
            dump(IterativeID_correlations, open(f'modeling_files/correlations/IterativeID_correlations_{_name}_{dataName}.json', 'w'))
            ASVset_abundances_df.index = [iterativeIDs_2.get(ASV, ASV).split('.')[0] for ASV in ASVset_abundances_df.index]
            ASVset_abundances_df.groupby(ASVset_abundances_df.index).sum()
            ASVset_correlations = {}
            for genus, sampleAbun in ASVset_abundances_df.iterrows():
                correlation, _p = _correlations(sampleAbun, Series(sample_values))
                if isna(correlation):
                    continue
                ASVset_correlations[genus] = {'correlation': correlation, 'p_value': _p}
            ASVset_correlations = dict(sorted(ASVset_correlations.items(), key=lambda item: item[1]['correlation'], reverse=True))
            dump(ASVset_correlations, open(f'modeling_files/correlations/ASVset_correlations_{_name}_{dataName}.json', 'w'))
    return ASVset_correlations, IterativeID_correlations, ommitted


@app.cell
def _(ASVset_correlations, IterativeID_correlations, display):
    display(IterativeID_correlations)
    print(len(IterativeID_correlations), len(ASVset_correlations))
    return


@app.cell
def _(DataFrame, Series, abs, display, load, multipletests, read_csv, sns):
    from matplotlib import colors, patches
    from pandas import set_option
    from glob import glob
    set_option('display.max_rows', None)
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    _sample_days = {'10AB': 0, 'I34': 149, 'J12': 161, 'K12': 182, 'L12': 189, 'M12': 210, 'N12': 238, 'P12': 283, 'Q': 300}
    order = ['all samples $H_2$ feed', 'all samples $H_2$ BT']
    _correlations, _pvals = ({}, {})
    for k in order:
        _correlations.setdefault(k, {})
        _pvals.setdefault(k, {})
    for cor in glob('modeling_files/correlations/ASVset_correlations_all_relevant_samples_H2*.json'):
        _name = cor.split('/')[-1].split('.')[0].split('_correlations_')[1].replace('_', ' ').replace('relevant ', '').replace('phase ', '').replace('H2', '$H_2$')
        content = load(open(cor, 'r'))
        _correlations[_name].update({k: v['correlation'] for k, v in content.items()})
        _pvals[_name].update({k: v['p_value'] for k, v in content.items()})
    q_vals = {}
    for _name, content in _pvals.items():
        pvals_list = list(content.values())
        q_vals[_name] = multipletests(pvals_list, alpha=0.05, method='fdr_bh')
        print(sum(q_vals[_name][0]))
    new_correlations = {_name: dict(inner) for _name, inner in _correlations.items()}
    new_pvals = {_name: dict(inner) for _name, inner in _pvals.items()}
    for _name, content in _pvals.items():
        for _i, (k, v) in enumerate(content.items()):
            if not q_vals[_name][0][_i]:
                new_correlations[_name].pop(k)
                new_pvals[_name].pop(k)
    display(new_correlations)
    display(new_pvals)
    df_2 = DataFrame(new_correlations).fillna(0)
    _pval_matrix = DataFrame(new_pvals)
    df_2.rename(columns={'all samples $H_2$ feed': '$H_2$ loading $\\left(\\frac{mol}{min}\\right)$', 'all samples $H_2$ BT': '$H_2$ breakthrough $\\left(\\frac{mol}{min}\\right)$'}, inplace=True)
    _pval_matrix.rename(columns={'all samples $H_2$ feed': '$H_2$ loading', 'all samples $H_2$ BT': '$H_2$ breakthrough'}, inplace=True)
    reduced = True
    if reduced:
        df_2 = df_2[(df_2 != float(0)).any(axis=1)]
    _total_df = read_csv('model_inputs/total.csv').set_index('seq')
    _orgs = {}
    for _i in df_2.index:
        _orgs.setdefault(_i.split('.')[0], []).append(_i)
    _taxonomies = {}
    level_2 = 'Genus'
    for _seq, _row in _total_df.iterrows():
        _day = _sample_days.get(_row['sample'])
        taxonomy_2 = []
        for l in reversed(_taxonomic_levels):
            _taxa = str(_row[l])
            IDs = _orgs.get(_taxa)
            if IDs is not None:
                break
        if _day is None or _row['rel_ab'] == 0 or IDs is None:
            continue
        taxonomy_string = '|'.join([str(_row[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index(level_2)])
        for _ID in IDs:
            _taxonomies.setdefault(_ID, taxonomy_string)
    _taxonomy_series = Series({_ID: _taxonomies.get(_ID, f"Unknown|{_ID.split('.')[0]}") for _ID in df_2.index})
    min_max = (round(df_2.min().min(), 1), round(df_2.max().max(), 1))
    _clusterMap = sns.clustermap(df_2, cmap='coolwarm_r', norm=colors.TwoSlopeNorm(vmin=min_max[0], vcenter=0, vmax=min_max[1]), col_cluster=False, clip_on=True, figsize=(20, 80), cbar_kws={'label': 'Correlation'})
    _labelsize = 70
    _clusterMap.ax_heatmap.set_xlabel('Operational metric', fontsize=_labelsize, labelpad=40)
    _clusterMap.ax_heatmap.set_ylabel('ASV', fontsize=_labelsize, labelpad=40)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_row_dendrogram.text(0.6, 0, 'Taxonomy Tree', transform=_clusterMap.ax_row_dendrogram.transAxes, fontsize=25, ha='center', va='top')
    _cbar = _clusterMap.ax_cbar
    _ticks = _cbar.get_yticks()[::2]
    _ticks[-1] = min_max[1]
    _ticks[0] = min_max[0]
    _ticks = [round(t, 1) for t in _ticks]
    _cbar.set_yticklabels(_ticks, fontsize=40)
    print(min_max, _ticks)
    _cbar.set_yticks(_ticks)
    _cbar.set_xlabel('Spearman $\\rho$', fontsize=60, labelpad=30)
    _clusterMap.ax_cbar.set_ylabel('')
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.55, _heatmap_pos.y0 + 0.45, 0.08, _heatmap_pos.height / 8])
    if reduced:
        _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 + 0.08, _heatmap_pos.y0 - 0.12, 0.1, _heatmap_pos.height / 7])
    _organisms_to_highlight = ['Methanobacterium', 'Methanosarcina', 'Methanobacteriaceae']
    iterativeID_levels_1 = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {k.split('.')[0]: v for k, v in iterativeID_levels_1.items()}
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    for _label in _ylabels:
        if any([x in _label.get_text() for x in _organisms_to_highlight]):
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
        _label.set_rotation(90)
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = list(range(_pval_matrix.shape[1]))
    pvals_reordered = _pval_matrix.iloc[_dendrogram_row, _dendrogram_col]
    ax = _clusterMap.ax_heatmap
    for _i in range(pvals_reordered.shape[0]):
        for _j in range(pvals_reordered.shape[1]):
            _p = pvals_reordered.iloc[_i, _j]
            _value = _clusterMap.data2d.iloc[_i, _j]
            if _value == 0:
                continue
            marker = ''
            _color = 'black' if abs(_value) < 0.7 else 'white'
            _clusterMap.ax_heatmap.text(_j + 0.5, _i + 0.5, f'{_value:.2f}{marker}', ha='center', va='center', color=_color, fontsize=60)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=50, rotation=80)
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=50, rotation=0)
    _clusterMap.figure.savefig(f"correlations_heatmap{('' if not reduced else '_reduced')}.png", bbox_inches='tight', dpi=300)
    return df_2, patches


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # All v All correlation matrix
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### Abundance
    """)
    return


@app.cell
def _(
    DataFrame,
    Patch,
    Series,
    df_2,
    display,
    json,
    load,
    nan,
    patches,
    read_csv,
    sns,
    spearmanr,
):
    from numpy import ones, triu, ones_like

    def corr_pvalues(df):
        n = df.shape[1]
        _pvals = DataFrame(ones((n, n)), index=df.columns, columns=df.columns)
        corrs = DataFrame(ones((n, n)), index=df.columns, columns=df.columns)
        for _i in range(n):
            for _j in range(_i + 1, n):
                c, _p = spearmanr(df.iloc[:, _i], df.iloc[:, _j])
                _pvals.iloc[_i, _j] = _p
                corrs.iloc[_i, _j] = c
                _pvals.iloc[_j, _i] = _p
                corrs.iloc[_j, _i] = c
        return (_pvals, corrs)
    _significantly_connected_organisms = [str(x) for x in load('FDR_passing_pairs.npy')]
    _significantly_connected_organisms.append('Methanobacteriaceae.1')
    print(_significantly_connected_organisms)
    _iterativeID_color_map = json.load(open(f'iterativeID_color_map.json', 'r'))
    iterativeIDs_3 = json.load(open('iterativeIDs.json', 'r'))
    _abundances = read_csv('abundances.csv', header=0).set_index('sample')
    display(_abundances)
    _abundances.drop([col for col in _abundances.columns if col not in _significantly_connected_organisms], axis=1, inplace=True)
    _abundances.drop([col for col in _abundances.columns if _abundances[col].max() < 0.005], axis=1, inplace=True)
    _total_captured = _abundances.sum(axis=1)
    print(_total_captured)
    display(_abundances.head())
    _pval_matrix, corr_matrix = corr_pvalues(_abundances)
    display(corr_matrix.head())
    display(_pval_matrix.head())
    taxonomy_3 = json.load(open('iterativeID_taxonomy.json', 'r'))
    _taxonomies = {}
    for col in _abundances.columns:
        _taxonomies[col] = '|'.join([v for k, v in taxonomy_3.get(col, 'Unknown').items() if k != 'Species' and v is not None])
    _taxonomy_series = Series({idx: _taxonomies.get(idx, f'Unknown|{idx}') for idx in _abundances.columns})
    print(f"Detected depth: {max((len(t.split('|')) for t in _taxonomy_series))}")
    print(f'df rows: {len(df_2.columns)}')
    print(f'taxonomy_series length: {len(_taxonomy_series)}')
    print(f'Sample entries:\n{_taxonomy_series.head()}')
    taxa_color_map_1 = json.load(open('iterativeID_color_map.json'))
    DEFAULT_COLOR = 'lightgray'
    bar_label = 'Phylum'
    row_colors = Series({idx: taxa_color_map_1.get(idx, DEFAULT_COLOR) for idx in corr_matrix.index}, name=bar_label)
    col_colors = Series({col: taxa_color_map_1.get(col, DEFAULT_COLOR) for col in corr_matrix.columns}, name=bar_label)
    _clusterMap = sns.clustermap(corr_matrix, row_colors=row_colors, col_colors=col_colors, cmap='coolwarm_r', center=0, figsize=(60, 70), dendrogram_ratio=(0.1, 0.2))
    _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    secondary_labels = [tick.get_text() for tick in _clusterMap.ax_heatmap.get_yticklabels()]
    _clusterMap.ax_heatmap.yaxis.set_ticks_position('left')
    _clusterMap.ax_heatmap.yaxis.set_label_position('left')
    _clusterMap.ax_heatmap.set_yticklabels(secondary_labels, rotation=0)
    hm_pos = _clusterMap.ax_heatmap.get_position()
    fig_w = _clusterMap.figure.get_figwidth()
    fig_h = _clusterMap.figure.get_figheight()
    strip_w = 0.015
    strip_h = 0.015
    _clusterMap.ax_row_colors.set_position([hm_pos.x0 - strip_w, hm_pos.y0, strip_w, hm_pos.height])
    _clusterMap.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0 - strip_h, hm_pos.width, strip_h])
    y_pad = strip_w * fig_w * 72 + 15
    _clusterMap.ax_heatmap.tick_params(axis='y', pad=y_pad)
    x_pad = strip_h * fig_h * 72 + 15
    _clusterMap.ax_heatmap.tick_params(axis='x', pad=x_pad)
    _cbar = _clusterMap.ax_cbar
    _ticks = _cbar.get_yticks()
    _ticks[-1] = corr_matrix.max().max()
    _ticks[0] = round(corr_matrix.min().min(), 2)
    _cbar.set_yticklabels(_ticks, fontsize=45)
    _cbar.set_yticks(_ticks)
    _cbar.set_xlabel('Spearman Correlation', fontsize=45, labelpad=30)
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.96, _heatmap_pos.y0 + 0.68, 0.06, _heatmap_pos.height / 5])
    _labelsize = 40
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=_labelsize, rotation=0)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=_labelsize, rotation=60, ha='right', rotation_mode='anchor')
    iterativeID_levels_2 = json.load(open('iterativeID_levels.json', 'r'))
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    _orgs = set()
    for _label in _ylabels:
        _text = _label.get_text()
        if iterativeID_levels_2.get(_text) == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels)
    _ID_levels = {k.split('.')[0]: v for k, v in iterativeID_levels_2.items()}
    _xlabels = _clusterMap.ax_heatmap.get_xticklabels()
    for _label in _xlabels:
        _text = _label.get_text()
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_xticklabels(_xlabels)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = _clusterMap.dendrogram_col.reordered_ind
    _one_triangle = True
    if _one_triangle:
        _df_reordered = corr_matrix.iloc[_dendrogram_row, _dendrogram_col]
        _mask = triu(ones_like(_df_reordered, dtype=bool), k=1)
        _mesh = _clusterMap.ax_heatmap.collections[0]
        _arr = _mesh.get_array().reshape(_df_reordered.shape)
        _arr[_mask] = nan
        _mesh.set_array(_arr.ravel())
        _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.2, _heatmap_pos.y0 + 0.4, 0.06, _heatmap_pos.height / 5])
    for _org in _orgs:
        _orgIx = corr_matrix.index.get_loc(_org)
        if _orgIx not in _dendrogram_row:
            continue
        _row_pos = _dendrogram_row.index(_orgIx)
        _rect = patches.Rectangle((0, _row_pos), len(corr_matrix.columns) if not _one_triangle else _row_pos + 1, 1, linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
        _colIx = corr_matrix.columns.get_loc(_org)
        _col_pos = _dendrogram_col.index(_colIx)
        _rect = patches.Rectangle((_col_pos, 0) if not _one_triangle else (_col_pos, len(corr_matrix.index)), 1, len(corr_matrix.index) if not _one_triangle else -(len(corr_matrix.index) - _col_pos), linewidth=6, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
    if isinstance(row_colors, Series):
        rc = row_colors
    else:
        rc = Series(list(row_colors), index=df_2.index)
    phylum_color = {}
    for idx in corr_matrix.index:
        parts = str(_taxonomies.get(idx, '')).split('|')
        if len(parts) < 2:
            continue
        _phylum = parts[1]
        if _phylum in ('None', '', 'Unknown'):
            continue
        _color = rc.get(idx)
        if _color is None:
            continue
        phylum_color.setdefault(_phylum, _color)
    print('phylum_color', phylum_color)
    archaea_markers = ('archaeo', 'halobacterota', 'methanobacteriota')

    def is_archaea(p):
        return any((m in _p.lower() for m in archaea_markers))
    archaea = sorted((_p for _p in phylum_color if is_archaea(_p)))
    bacteria = sorted((_p for _p in phylum_color if not is_archaea(_p)))
    handles = []
    if archaea:
        handles.append(Patch(color='none', label='$\\bf{Archaea}$'))
        handles = handles + [Patch(facecolor=phylum_color[_p], label=_p) for _p in archaea]
    if bacteria:
        handles.append(Patch(color='none', label='$\\bf{Bacteria}$'))
        handles = handles + [Patch(facecolor=phylum_color[_p], label=_p) for _p in bacteria]
    _clusterMap.figure.legend(handles=handles, title='Phylum', title_fontsize=50, fontsize=40, loc='upper right', bbox_to_anchor=(0.6, 0.7), frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6)
    _clusterMap.ax_heatmap.set_xlabel('Member ASVs', fontsize=50)
    _clusterMap.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)
    _clusterMap.figure.savefig(f"abundance_correlatons{('_one_triangle' if _one_triangle else '')}.png", dpi=300, bbox_inches='tight')
    return corr_matrix, ones, ones_like, taxa_color_map_1, triu


@app.cell
def _(KMeans, best_k, corr_matrix, pd, sig_mask):
    # sig_mask is a boolean DataFrame (True = p < threshold)
    corr_masked = corr_matrix.where(sig_mask, 0)
    km = KMeans(n_clusters=best_k, n_init=20, random_state=42)
    cluster_labels = pd.Series(km.fit_predict(corr_masked), index=corr_matrix.index)
    cluster_labels
    return


@app.cell
def _(
    DataFrame,
    abs,
    display,
    inf,
    load,
    nan,
    ones,
    ones_like,
    patches,
    sns,
    spearmanr,
    triu,
):
    days = {'I34': 149, 'J12': 161, 'K12': 182, 'L12': 189, 'M12': 210, 'N12': 238, 'P12': 283, 'Q': 300}

    def corr_pvalues_1(df):
        n = df.shape[1]
        _pvals = DataFrame(ones((n, n)), index=df.columns, columns=df.columns)
        for _i in range(n):
            for _j in range(_i + 1, n):
                _, _p = spearmanr(df.iloc[:, _i], df.iloc[:, _j])
                _pvals.iloc[_i, _j] = _p
                _pvals.iloc[_j, _i] = _p
        return _pvals
    _nonzero_per_day = load(open('model_inputs/nonzero_per_day.json', 'r'))
    _top_per_day = {}
    _all_orgs = {}
    _topNum = 10
    for _day, _org_dict in _nonzero_per_day.items():
        _orgs = dict(list(_org_dict.items())[:_topNum])
        _top_per_day[_day] = _orgs
        _all_orgs.update(_orgs)
    df_3 = DataFrame(_top_per_day).T
    df_3 = df_3.astype(float).replace([inf, -inf], nan)
    df_3 = df_3.loc[:, df_3.notna().sum() >= 5]
    display(df_3)
    corr_matrix_1 = df_3.corr('spearman').dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    display(corr_matrix_1)
    _pval_matrix = corr_pvalues_1(df_3).dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    _pval_matrix = _pval_matrix[corr_matrix_1.columns]
    _pval_matrix = _pval_matrix.loc[corr_matrix_1.index]
    display(_pval_matrix)
    _clusterMap = sns.clustermap(corr_matrix_1, cmap='coolwarm_r', center=0, figsize=(60, 60), col_cluster=True, row_cluster=True, dendrogram_ratio=(0.1, 0.2))
    _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
    _dendrogram_row = _clusterMap.dendrogram_row.reordered_ind
    _dendrogram_col = _clusterMap.dendrogram_col.reordered_ind
    _one_triangle = True
    if _one_triangle:
        _df_reordered = corr_matrix_1.iloc[_dendrogram_row, _dendrogram_col]
        _mask = triu(ones_like(_df_reordered, dtype=bool), k=1)
        _mesh = _clusterMap.ax_heatmap.collections[0]
        _arr = _mesh.get_array().reshape(_df_reordered.shape)
        _arr[_mask] = nan
        _mesh.set_array(_arr.ravel())
    _organisms_to_highlight = ['Methanobacterium.2', 'Methanobacteriaceae.1', 'Methanobacterium.1']
    for _org in _organisms_to_highlight:
        _orgIx = corr_matrix_1.index.get_loc(_org)
        if _orgIx not in _dendrogram_row:
            continue
        _row_pos = _dendrogram_row.index(_orgIx)
        _rect = patches.Rectangle((0, _row_pos), len(corr_matrix_1.columns) if not _one_triangle else _row_pos + 1, 1, linewidth=16, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
        _colIx = corr_matrix_1.columns.get_loc(_org)
        _col_pos = _dendrogram_col.index(_colIx)
        _rect = patches.Rectangle((_col_pos, 0) if not _one_triangle else (_col_pos, len(corr_matrix_1.index)), 1, len(corr_matrix_1.index) if not _one_triangle else -(len(corr_matrix_1.index) - _col_pos), linewidth=16, edgecolor='black', facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_rect)
    _cbar = _clusterMap.ax_cbar
    _ticks = _cbar.get_yticks()
    _ticks[-1] = round(corr_matrix_1.max().max())
    _ticks[0] = round(corr_matrix_1.min().min(), 2)
    _cbar.set_yticklabels(_ticks, fontsize=50)
    _cbar.set_yticks(_ticks)
    _cbar.set_ylabel('Spearman $\\rho$', fontsize=60, labelpad=15, rotation=90)
    _labelsize = 60
    _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=_labelsize, rotation=0)
    _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=_labelsize, rotation=60, ha='right', rotation_mode='anchor')
    for _i in range(_clusterMap.data2d.shape[0]):
        for _j in range(_clusterMap.data2d.shape[1]):
            if _mask[_i, _j] == 1:
                continue
            _value = _clusterMap.data2d.iloc[_i, _j]
            if _value == 0:
                continue
            _color = 'black' if abs(_value) < 0.7 else 'white'
            _clusterMap.ax_heatmap.text(_j + 0.5, _i + 0.5, f'{_value:.2f}', ha='center', va='center', color=_color, fontsize=60)
    iterativeID_levels_3 = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {k.split('.')[0]: v for k, v in iterativeID_levels_3.items()}
    _clusterMap.ax_heatmap.yaxis.set_ticks_position('left')
    _clusterMap.ax_heatmap.yaxis.set_label_position('left')
    _ylabels = _clusterMap.ax_heatmap.get_yticklabels()
    for _label in _ylabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_yticklabels(_ylabels, rotation=0)
    _xlabels = _clusterMap.ax_heatmap.get_xticklabels()
    for _label in _xlabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_xticklabels(_xlabels)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 1.01, _heatmap_pos.y0 + 0.5, 0.06, _heatmap_pos.height / 5])
    if _one_triangle:
        _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.3, _heatmap_pos.y0 + 0.4, 0.06, _heatmap_pos.height / 5])
    _clusterMap.ax_heatmap.set_xlabel('ASVs', fontsize=100, labelpad=20)
    _clusterMap.ax_heatmap.set_ylabel('ASVs', fontsize=100, labelpad=20)
    _clusterMap.figure.savefig(f"abundance_heatmaps/Top_{_topNum}_ASVs_abundance_correlation{('_one_triangle' if _one_triangle else '')}.png", bbox_inches='tight', dpi=300)
    return (corr_pvalues_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # co-occurrence figure
    """)
    return


@app.cell
def _(
    archaea_phyla,
    bacteria_phyla,
    combinations,
    corr_pvalues_1,
    defaultdict,
    display,
    json,
    load,
    mcolors,
    mpatches,
    np,
    order_color_map,
    pair_data,
    plt,
    read_csv,
    reject,
    taxa_color_map_1,
):
    import networkx as nx
    import matplotlib.cm as cm
    import math
    iterativeIDs_4 = load(open('iterativeIDs.json'))
    _iterativeID_taxonomy = load(open('iterativeID_taxonomy.json'))
    level_3 = 'Phylum'
    iterativeID_level = {_ID: content.get(level_3, 'Unknown') for _ID, content in _iterativeID_taxonomy.items()}
    genera_level = False
    if genera_level:
        iterativeID_level = iterativeID_level = {_ID.split('.')[0]: k for _ID, k in iterativeID_level.items()}
    _significantly_connected_organisms = [str(x) for x in np.load('FDR_passing_pairs.npy')]
    _significantly_connected_organisms.append('Methanobacteriaceae.1')
    print(_significantly_connected_organisms)
    _iterativeID_color_map = json.load(open(f'iterativeID_color_map.json', 'r'))
    iterativeIDs_4 = json.load(open('iterativeIDs.json', 'r'))
    _abundances = read_csv('abundances.csv', header=0).set_index('sample')
    display(_abundances)
    _abundances.drop([col for col in _abundances.columns if col not in _significantly_connected_organisms], axis=1, inplace=True)
    pvals_matrix, corr_matrix_2 = corr_pvalues_1(_abundances)
    _total_captured = _abundances.sum(axis=1)
    print(_total_captured)
    display(_abundances.head())
    _mean_rel_abund = _abundances.mean(axis=0)
    _presence = (_abundances > 0).astype(int)
    _cooccurrence = defaultdict(int)
    for _sample in _presence.itertuples(index=False):
        _present = [col for col, val in zip(_presence.columns, _sample) if val]
        for _pair in combinations(sorted(_present), 2):
            _cooccurrence[_pair] = _cooccurrence[_pair] + 1
    G = nx.Graph()
    for (_a, _b), _count in _cooccurrence.items():
        _rho = corr_matrix_2.loc[_a, _b]
        _p_value = pvals_matrix.loc[_a, _b]
        if np.isnan(_rho):
            continue
        G.add_edge(_a, _b, weight=np.abs(_rho), rho=_rho, pvalue=_p_value, cooccurrence=_count)
    print(f'Tests run: {len(pair_data)}')
    print(f'FDR-significant pairs (q < 0.05): {int(reject.sum())}')
    print(f'Edges after FDR correction: {G.number_of_edges()}')
    print(f'Nodes in graph: {G.number_of_nodes()}')
    print('nodes:', G.number_of_nodes(), 'edges:', G.number_of_edges())
    print('max cooccurrence:', max(_cooccurrence.values()) if _cooccurrence else 0)

    def header_patch(title):
        return mpatches.Patch(color='none', label=f'$\\bf{{{title}}}$')

    def render_network(G_sub, title_slug, focus_nodes=None):
        if G_sub.number_of_edges() == 0:
            print(f'[{title_slug}] no edges — skipping')
            return
        pos = nx.spring_layout(G_sub, seed=42, iterations=100, k=0.3)
        edges = G_sub.edges(data=True)
        rho_values = [d['rho'] for _, _, d in edges]
        edge_widths = [3 * d['weight'] for _, _, d in edges]
        norm = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        cmap = cm.RdBu
        edge_colors = [cmap(norm(r)) for r in rho_values]
        width, height = (40, 30)
        fig, ax = plt.subplots(figsize=(width, height))
        scale = 5000 * (width / 10) * 2
        compressor = np.sqrt
        node_sizes = [scale * compressor(_mean_rel_abund.get(n, 0)) for n in G_sub.nodes()]
        node_colors = [taxa_color_map_1.get(iterativeID_level.get(n, 'Unknown'), 'lightgray') for n in G_sub.nodes()]
        print('unique node colors:', len(set(node_colors)))
        edgecolors = ['black' if focus_nodes and n in focus_nodes else 'none' for n in G_sub.nodes()]
        linewidths = [3 if focus_nodes and n in focus_nodes else 0 for n in G_sub.nodes()]
        nx.draw_networkx_nodes(G_sub, pos, node_size=node_sizes, node_color=node_colors, edgecolors=edgecolors, linewidths=linewidths, alpha=0.9, ax=ax)
        nx.draw_networkx_edges(G_sub, pos, width=edge_widths, edge_color=edge_colors, alpha=0.85, ax=ax)
        for n, (x, y) in pos.items():
            abund = _mean_rel_abund.get(n, 0)
            font_size = 6 + 14 * compressor(abund) / compressor(_mean_rel_abund.max()) * (width / 10)
            font_size = max(5, min(font_size, 48))
            ax.text(x, y, str(n), fontsize=font_size, color='black', fontweight='bold', ha='center', va='center')
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        _cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        _cbar.set_label('Spearman ρ', fontsize=10 * (width / 10))
        _cbar.ax.tick_params(labelsize=8 * (width / 10), length=8, width=2)
        legend_entries = [(0.002, '0.2%'), (0.02, '2%'), (0.2, '20%')]
        sizes_pt2 = [scale * compressor(_a) for _a, _ in legend_entries]
        diameters_pt = [2 * np.sqrt(s / np.pi) for s in sizes_pt2]
        breather_pt = 8
        offsets_pt = [0.0]
        for _i in range(1, len(legend_entries)):
            offsets_pt.append(offsets_pt[-1] + (diameters_pt[_i - 1] + diameters_pt[_i]) / 2 + breather_pt)
        fig = ax.figure
        fig_h_pts = fig.get_figheight() * 72
        top_y = 0.9
        circle_x = 0.85
        label_x = 0.88
        ax.text(circle_x, top_y + 0.025, 'Mean rel. abundance', transform=fig.transFigure, va='bottom', fontweight='bold', fontsize=6 * (width / 10), clip_on=False)
        for (abund, _label), s, off in zip(legend_entries, sizes_pt2, offsets_pt):
            y = top_y - off / fig_h_pts
            ax.scatter([circle_x], [y], s=s, color='slategray', alpha=0.9, transform=fig.transFigure, clip_on=False)
            ax.text(label_x, y, _label, transform=fig.transFigure, va='center', fontsize=5 * (width / 10), clip_on=False)
        archaea_patches = [mpatches.Patch(color=order_color_map[_p], label=_p) for _p in archaea_phyla if _p in order_color_map]
        bacteria_patches = [mpatches.Patch(color=order_color_map[_p], label=_p) for _p in bacteria_phyla if _p in order_color_map]
        legend_handles = [header_patch('Archaea')] + archaea_patches + [header_patch('Bacteria')] + bacteria_patches
        ax.legend(handles=legend_handles, title='Taxonomic ' + level_3, title_fontsize=8 * (width / 10), loc='lower left', bbox_to_anchor=(-0.24, 0.1), fontsize=7 * (width / 10), frameon=True)
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(f'cooccurrence_network_{title_slug}.png', dpi=300, bbox_inches='tight')
        plt.show()
        plt.close(fig)
    selections = {'GAOs': ['Methanobacteriaceae.1', 'Methanosarcina.3', ...], 'PAOs': ['Anaerolinea.13', ...]}
    for _name, focus in selections.items():
        focus_in_G = [n for n in focus if n in G]
        if not focus_in_G:
            print(f'[{_name}] none of the focus nodes are in G — skipping')
            continue
        keep = set(focus_in_G)
        for n in focus_in_G:
            keep.update(G.neighbors(n))
        G_sub = G.subgraph(keep).copy()
        render_network(G_sub, _name, focus_nodes=set(focus_in_G))
    return


if __name__ == "__main__":
    app.run()
