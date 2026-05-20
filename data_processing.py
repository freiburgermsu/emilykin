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
def _(mo):
    from json import load, dump
    from pandas import read_csv
    _total_df = read_csv('table_rel_export.csv').set_index('seq')
    _total_df = _total_df.fillna('')
    mo.output.append(_total_df)
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
def _(dump, mo, taxonomy):
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    iterativeIDs, _iterativeID_levels, iterativeID_phylums, taxonIDs, iterativeID_taxonomy = ({}, {}, {}, {}, {})
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
        _iterativeID_levels[taxon] = best_level
        iterativeID_phylums[taxon] = _taxa['Phylum']
        iterativeID_taxonomy[taxon] = taxonomy.get(_seq, 'Unknown')
    dump(iterativeIDs, open('iterativeIDs.json', 'w'))
    dump(_iterativeID_levels, open('iterativeID_levels.json', 'w'))
    dump(iterativeID_phylums, open('iterativeID_phylums.json', 'w'))
    dump(taxonIDs, open('taxonIDs.json', 'w'))
    dump(iterativeID_taxonomy, open('iterativeID_taxonomy.json', 'w'))
    print(len(iterativeIDs), 'unique orgs')
    print(len(taxonIDs), 'unique taxa')
    mo.output.append(list(iterativeIDs.items())[-10:])
    mo.output.append(list(_iterativeID_levels.items())[-10:])
    mo.output.append(list(iterativeID_phylums.items())[-10:])
    return iterativeIDs, level


@app.cell
def _(iterativeIDs, load, mo):
    from pandas import DataFrame
    _ab = load(open('relative_abundance.json', 'r'))
    ab_df = DataFrame(_ab) / 100
    _seqIDs = {_seq: _ID for _ID, _seq in iterativeIDs.items()}
    ab_df.columns = [_seqIDs[col] for col in ab_df.columns]
    print(ab_df.shape, min(ab_df.min()), max(ab_df.max()))
    not_normalized = [x for x, val in ab_df.T.sum().items() if round(val, 1) != 1]
    mo.output.append(len(not_normalized)); mo.output.append(not_normalized)
    _zero_level = 1e-05
    ab_df = ab_df.loc[:, (ab_df.fillna(0) > _zero_level).sum() >= 6]
    mo.output.append(ab_df.shape)
    mo.output.append(ab_df.head())
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
def _(defaultdict, np):
    def taxonomy_linkage(taxonomy_series):
        """
        Build a linkage matrix that EXACTLY follows taxonomy hierarchy.
        Auto-detects depth from the taxonomy strings.
        """
        parsed = taxonomy_series.fillna('unknown').astype(str).str.split('[;|,]', regex=True).apply(lambda lst: [x.strip() for x in lst or [] if x.strip()])
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
        Z = np.array(linkage_rows, dtype=float)
        if len(Z) > 0:
            for _i in range(1, len(Z)):
                if Z[_i, 2] < Z[_i - 1, 2]:
                    Z[_i, 2] = Z[_i - 1, 2]
        return Z
    GAOs_PAOs = {'PAOs': ['Ca_Accumulibacter', 'Tetrasphaera', 'Dechloromonas', 'Microlunatus', 'Azonexus', 'Ca_Phosphoribacter'], 'GAOs': ['Ca_Competibacter', 'Defluviicoccus', 'Propionivibrio', 'Ca_Contendobacter'], 'Putative PAOs': ['Ca_Obscuribacter', 'Thauera', 'Zoogloea', 'Paracoccus'], 'Putative GAOs': ['Micropruina', 'Amaricoccus', 'Ca_Glycocaulis', 'Thauera'], 'Other PHA storing potential+ function': ['Pseudomonas', 'Bacillus', 'Acinetobacter', 'Rhodocyclaceae']}
    return GAOs_PAOs, taxonomy_linkage


@app.cell
def _(dump, level, load):
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from collections import defaultdict
    from itertools import combinations
    import numpy as np
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import multipletests

    def rgba_to_hex(rgba):
        return '#{:02x}{:02x}{:02x}'.format(*(int(c * 255) for c in rgba[:3]))
    iterativeID_taxonomy_1 = load(open('iterativeID_taxonomy.json', 'r'))
    iterativeID_phylums_1 = load(open('iterativeID_phylums.json', 'r'))
    none_keys = [_k for _k, _v in iterativeID_phylums_1.items() if _v is None]
    print('None taxonomy entries:', none_keys)
    archaea_phyla = sorted({_v for _k, _v in iterativeID_phylums_1.items() if _v is not None and ('archaeo' in _v.lower() or any((_a in _v.lower() for _a in ['candidatus thermoplasmatota', 'halobacterota', 'methanobacteriota'])))})
    bacteria_phyla = sorted({_v for _k, _v in iterativeID_phylums_1.items() if _v is not None and _v not in archaea_phyla})
    taxa_color_map = {}
    for _i, _p in enumerate(archaea_phyla):
        taxa_color_map[_p] = rgba_to_hex(plt.cm.turbo(_i / max(len(archaea_phyla), 1) * 0.15))
    for _i, _p in enumerate(bacteria_phyla):
        taxa_color_map[_p] = rgba_to_hex(plt.cm.turbo(0.2 + _i / max(len(bacteria_phyla), 1) * 0.8))
    iterativeID_color_map = {_ID: taxa_color_map[phylum] for _ID, phylum in iterativeID_phylums_1.items() if phylum}

    def expand_phylum_to_classes(phylum, cmap, taxonomy, iterativeID_color_map, taxa_color_map, lo=0.6, hi=0.95, base_t=0.85):
        classes = sorted({t['Class'] for t in taxonomy.values() if t.get('Phylum') == phylum and t.get('Class')})
        n = max(len(classes), 1)
        class_colors = {c: rgba_to_hex(cmap(lo + (hi - lo) * _i / max(n - 1, 1))) for _i, c in enumerate(classes)}
        base_color = rgba_to_hex(cmap(base_t))
        for org_id, _taxa in taxonomy.items():
            if _taxa.get('Phylum') == phylum and _taxa.get('Class') in class_colors:
                iterativeID_color_map[org_id] = class_colors[_taxa['Class']]
        taxa_color_map[phylum] = base_color
        return (class_colors, base_color)
    proteo_class_color, _proteo_base = expand_phylum_to_classes('Proteobacteria', plt.cm.Purples, iterativeID_taxonomy_1, iterativeID_color_map, taxa_color_map)
    dump(taxa_color_map, open(f'{level}_color_map.json', 'w'))
    dump(iterativeID_color_map, open('iterativeID_color_map.json', 'w'))
    dump(proteo_class_color, open('proteo_class_color.json', 'w'))
    dump({'Proteobacteria': _proteo_base}, open('phylum_base_overrides.json', 'w'))
    phylum_colors = load(open('phylum_colors.json', 'r'))
    DEFAULT = 'lightgray'
    genera_color_map = {ID.split('.')[0]: v for ID, v in iterativeID_color_map.items()}
    return (
        DEFAULT,
        archaea_phyla,
        bacteria_phyla,
        combinations,
        defaultdict,
        genera_color_map,
        iterativeID_color_map,
        iterativeID_phylums_1,
        iterativeID_taxonomy_1,
        mcolors,
        multipletests,
        np,
        plt,
        spearmanr,
    )


@app.cell
def _(
    DEFAULT,
    archaea_phyla,
    bacteria_phyla,
    colorsys,
    combinations,
    defaultdict,
    iterativeID_color_map,
    iterativeID_phylums_1,
    iterativeIDs,
    mcolors,
    multipletests,
    np,
    phylum_color_map,
    read_csv,
    spearmanr,
):
    import matplotlib.patches as mpatches
    df = read_csv('abundances.csv').set_index('sample')
    df = df.loc[:, (df.fillna(0) > 0).sum() >= 6]
    _seqIDs = {_k: _v for _v, _k in iterativeIDs.items()}
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
    member_colors = {iterativeID_phylums_1.get(n, n): iterativeID_color_map.get(n, DEFAULT) for n in mems}
    archaea_phyla_1 = {_k for _k in archaea_phyla if _k in member_colors.keys()}
    bacteria_phyla_1 = {_k for _k in bacteria_phyla if _k in member_colors.keys()}
    archaea_patches = [mpatches.Patch(color=member_colors[_p], label=_p) for _p in archaea_phyla_1]
    bacteria_patches = [mpatches.Patch(color=member_colors[_p], label=_p) for _p in bacteria_phyla_1]

    def _lighten(c, factor):
        h, l, s = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
        return colorsys.hls_to_rgb(h, max(0.0, min(1.0, l * factor)), s)
    _proteo_base = phylum_color_map.get('Proteobacteria', 'tab:purple')
    return mpatches, pair_data, reject


@app.cell
def _(
    DataFrame,
    GAOs_PAOs,
    genera_color_map,
    mo,
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
    iterativeID_color_map_1 = json.load(open('iterativeID_color_map.json'))
    phylum_color_map = json.load(open('Phylum_color_map.json'))
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    abundances = read_csv('abundances.csv').set_index('sample')
    abundances.index = [_sample_days[col] for col in abundances.index]
    abundances = abundances.loc[sorted(abundances.index, key=int)]
    intervals = [0] + list(logspace(-2, -0.3, 5))
    print(intervals)

    def create_heatmap(df, taxonomies, title, inlayed_data=None, linkage=None, mode='abundance'):
        if mode == 'log2fc':
            new_cmap = LinearSegmentedColormap.from_list('FCMap', [(0.0, '#2166ac'), (0.5, 'white'), (1.0, '#b2182b')])
            new_cmap.set_bad('lightgray')
            _abs_max = float(np.nanmax(np.abs(df.values))) if df.size else 1.0
            if not np.isfinite(_abs_max) or _abs_max == 0:
                _abs_max = 1.0
            vmin = -_abs_max
            vmax = _abs_max
            vcenter = 0.0
            norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        else:
            new_cmap = LinearSegmentedColormap.from_list('NewMap', [(0.0, 'aliceblue'), (0.25, 'lightblue'), (1.0, 'navy')])
            new_cmap.set_bad('aliceblue')
            vmin = df.min().min()
            vmax = df.max().max()
            vcenter = log10(0.1)
            vcenter = min(max(vcenter, vmin + 1e-06), vmax - 1e-06)
            norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        proteo_class_color = json.load(open('proteo_class_color.json'))
        _proteo_base = json.load(open('phylum_base_overrides.json'))['Proteobacteria']

        def lookup(idx):
            if idx in iterativeID_color_map_1:
                return iterativeID_color_map_1[idx]
            if idx in genera_color_map:
                return genera_color_map[idx]
            return _DEFAULT_COLOR
        row_colors = Series({idx: lookup(idx) for idx in df.index}, name='Phylum')
        if mode == 'log2fc':
            _figsize = (max(14, 1.5 * df.shape[1] + 8), 20)
            _dendro_ratio = (0.08, 0.15)
        else:
            _figsize = (50, 20)
            _dendro_ratio = (0.2, 0.15)
        _clusterMap = sns.clustermap(df, row_colors=row_colors, cmap=new_cmap, norm=norm, clip_on=True, col_cluster=False, row_cluster=True, figsize=_figsize, row_linkage=taxonomy_linkage(taxonomies), dendrogram_ratio=_dendro_ratio)
        for tick in _clusterMap.ax_row_colors.get_xticklabels():
            tick.set_fontsize(22)
        _cbar = _clusterMap.ax_cbar
        if mode == 'log2fc':
            _fc_ticks = np.linspace(vmin, vmax, 5)
            _cbar.set_yticks(_fc_ticks)
            _cbar.set_yticklabels([f'{t:+.1f}' for t in _fc_ticks], fontsize=22)
            _cbar.set_xlabel('log$_2$ fold change', fontsize=18, labelpad=10)
        else:
            log_ticks = np.array([vmin, log10(0.01), vcenter, log10(0.2), vmax])
            _cbar.set_yticks(log_ticks)
            original_ticks = 10 ** log_ticks
            original_labels = [f'{x * 100:.1e}'.replace('e-0', 'E-') if x * 100 < 0.1 else round(x * 100) for x in original_ticks]
            original_labels[0] = '0'
            _cbar.set_yticklabels(original_labels, fontsize=22)
            _cbar.set_yticks(log_ticks)
            _cbar.set_xlabel('Rel. Abundance %', fontsize=18, labelpad=10)
        # Shannon-diversity inlay disabled per request.
        # if inlayed_data:
        #     top_ax = _clusterMap.ax_col_dendrogram
        #     top_ax.clear()
        #     top_ax.plot(inlayed_data.keys(), inlayed_data.values(), color='black', marker='o', linewidth=3)
        #     top_ax.set_ylabel('Shannon Diversity', fontsize=20)
        #     top_ax.yaxis.tick_right()
        #     top_ax.yaxis.set_label_position('right')
        #     top_ax.grid(True, axis='y', linestyle='--', alpha=0.7)
        #     top_ax.set_xticks([])
        #     top_ax.tick_params(axis='y', labelsize=16)
        _clusterMap.ax_col_dendrogram.set_visible(False)
        _clusterMap.figure.subplots_adjust(bottom=0.15, top=0.95)
        _clusterMap.ax_heatmap.set_yticklabels(_clusterMap.ax_heatmap.get_yticklabels(), fontsize=24, rotation=0)
        _clusterMap.ax_heatmap.set_xticklabels(_clusterMap.ax_heatmap.get_xticklabels(), fontsize=24, rotation=70, ha='right')
        _xlabel = 'Phase comparison (log$_2$ fold change)' if mode == 'log2fc' else 'Days of operation'
        _clusterMap.ax_heatmap.set_xlabel(_xlabel, fontsize=32, labelpad=20)
        _clusterMap.ax_row_dendrogram.xaxis.set_visible(False)
        _clusterMap.ax_row_dendrogram.text(0.65, -0.02, 'Taxonomical tree', fontsize=24, ha='center', transform=_clusterMap.ax_row_dendrogram.transAxes)
        if mode == 'log2fc':
            # Layout: [phylum colors] [heatmap] [row labels] [compressed dendrogram]
            colors_left = 0.18
            colors_w = 0.018
            heatmap_left = colors_left + colors_w + 0.005
            heatmap_right = 0.58
            dendro_w = 0.06
            dendro_left = heatmap_right + 0.20
            dendro_right = dendro_left + dendro_w
            hm_pos = _clusterMap.ax_heatmap.get_position()
            dend_pos = _clusterMap.ax_row_dendrogram.get_position()
            _clusterMap.ax_heatmap.set_position([heatmap_left, hm_pos.y0, heatmap_right - heatmap_left, hm_pos.height])
            _clusterMap.ax_heatmap.yaxis.tick_right()
            _clusterMap.ax_heatmap.yaxis.set_label_position('right')
            _clusterMap.ax_heatmap.tick_params(axis='y', pad=4)
            # Tentative dendrogram placement; final x is set after labels render.
            _clusterMap.ax_row_dendrogram.set_position([heatmap_right + 0.20, dend_pos.y0, dendro_w, dend_pos.height])
            # Flip so branches open leftward, toward the heatmap.
            _clusterMap.ax_row_dendrogram.invert_xaxis()
            col_dend_pos = _clusterMap.ax_col_dendrogram.get_position()
            _cbar_w = 0.035
            _cbar_x = colors_left - _cbar_w - 0.04
            _clusterMap.ax_cbar.set_position([_cbar_x, col_dend_pos.y0, _cbar_w, col_dend_pos.height])
        else:
            _clusterMap.ax_heatmap.yaxis.tick_left()
            _clusterMap.ax_heatmap.yaxis.set_label_position('left')
            label_right = 0.22
            heatmap_right = 0.68
            dendro_left = 0.71
            dendro_right = 0.85
            hm_pos = _clusterMap.ax_heatmap.get_position()
            dend_pos = _clusterMap.ax_row_dendrogram.get_position()
            _clusterMap.ax_heatmap.set_position([label_right, hm_pos.y0, heatmap_right - label_right, hm_pos.height])
            _clusterMap.ax_row_dendrogram.set_position([dendro_left, dend_pos.y0, dendro_right - dendro_left, dend_pos.height])
            _clusterMap.ax_row_dendrogram.invert_xaxis()
            gap = 0.03
            col_dend_pos = _clusterMap.ax_col_dendrogram.get_position()
            _clusterMap.ax_col_dendrogram.set_position([label_right, col_dend_pos.y0 + gap, heatmap_right - label_right, col_dend_pos.height])
            _shannon_pos = _clusterMap.ax_col_dendrogram.get_position()
            _cbar_w = 0.05
            _cbar_x = _shannon_pos.x0 - _cbar_w - 0.02
            _clusterMap.ax_cbar.set_position([_cbar_x, _shannon_pos.y0, _cbar_w, _shannon_pos.height])
        _iterativeID_levels = json.load(open('iterativeID_levels.json', 'r'))
        _ID_levels = {}
        for _k, _v in _iterativeID_levels.items():
            _ID_levels[_k] = _v
            _ID_levels.setdefault(_k.split('.')[0], _v)
        inverted_GAOs_PAOs = {_v: _k for _k, vs in GAOs_PAOs.items() for _v in vs}
        for _label in _clusterMap.ax_heatmap.get_yticklabels():
            _text = _label.get_text()
            for _org, val in inverted_GAOs_PAOs.items():
                if _org not in _text:
                    continue
                _label.set_fontweight('bold')
                if 'GAOs' in val:
                    _label.set_color('green')
                    if 'Putative' in val:
                        _label.set_color('mediumseagreen')
                elif 'PAOs' in val:
                    _label.set_color('blue')
                    if 'Putative' in val:
                        _label.set_color('cornflowerblue')
                elif 'PHA' in val:
                    _label.set_color('red')
            if _ID_levels.get(_text) == 'Genus':
                _label.set_fontstyle('italic')
        hm_pos = _clusterMap.ax_heatmap.get_position()
        fig_w = _clusterMap.figure.get_figwidth()
        fig_h = _clusterMap.figure.get_figheight()
        strip_w = 0.015
        strip_h = 0.015
        if mode == 'log2fc':
            _clusterMap.ax_row_colors.set_position([hm_pos.x0 - strip_w - 0.005, hm_pos.y0, strip_w, hm_pos.height])
        else:
            _clusterMap.ax_row_colors.set_position([hm_pos.x1 + 0.005, hm_pos.y0, strip_w, hm_pos.height])
            _clusterMap.ax_heatmap.tick_params(axis='y', pad=4)
        if isinstance(row_colors, Series):
            rc = row_colors
        else:
            rc = Series(list(row_colors), index=df.index)
        phylum_color = {}
        for idx in df.index:
            parts = str(taxonomies.get(idx, '')).split('|')
            if len(parts) < 2:
                continue
            phylum = parts[1]
            if phylum in ('None', '', 'Unknown'):
                continue
            if phylum == 'Proteobacteria':
                phylum_color['Proteobacteria'] = _proteo_base
                continue
            _color = rc.get(idx)
            if _color is None:
                continue
            phylum_color.setdefault(phylum, _color)
        print('phylum_color', phylum_color)
        archaea_markers = ('archaeo', 'halobacterota', 'methanobacteriota')

        def is_archaea(p):
            return any((m in p.lower() for m in archaea_markers))
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
        if mode == 'log2fc':
            _ncol = int(np.ceil(len(handles) / 2))
            _clusterMap.figure.legend(handles=handles, title='Phylum', title_fontsize=20, fontsize=16, loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=_ncol, frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6, columnspacing=1.2)
        else:
            _clusterMap.figure.legend(handles=handles, title='Phylum', title_fontsize=20, fontsize=18, loc='upper right', bbox_to_anchor=(dendro_right, 1.0), frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6)
        for spine in _clusterMap.ax_heatmap.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor('black')
            spine.set_linewidth(1)
        try:
            _phase_day_ranges = json.load(open('phase_day_ranges.json', 'r'))
        except FileNotFoundError:
            _phase_day_ranges = {}
        if _phase_day_ranges and mode != 'log2fc':
            import matplotlib.lines as _mlines
            _days_int = []
            for _c in df.columns:
                try:
                    _days_int.append(int(_c))
                except (ValueError, TypeError):
                    _days_int = None
                    break
            if _days_int:
                _n_cols = len(_days_int)

                def _day_to_heatmap_x(_d):
                    if _d <= _days_int[0]:
                        return 0.5
                    if _d >= _days_int[-1]:
                        return _n_cols - 0.5
                    for _i in range(_n_cols - 1):
                        if _days_int[_i] <= _d <= _days_int[_i + 1]:
                            if _days_int[_i + 1] == _days_int[_i]:
                                return _i + 0.5
                            _frac = (_d - _days_int[_i]) / (_days_int[_i + 1] - _days_int[_i])
                            return _i + 0.5 + _frac
                    return None
                _fig = _clusterMap.figure
                _hm_pos = _clusterMap.ax_heatmap.get_position()
                _top_pos = _clusterMap.ax_col_dendrogram.get_position()
                _y_bot = _hm_pos.y0
                _y_top = _top_pos.y1
                _fig.canvas.draw()
                _lower_bounds = sorted({_span['start'] for _span in _phase_day_ranges.values()})
                _lower_bounds = _lower_bounds[1:]
                for _day in _lower_bounds:
                    _x_data = _day_to_heatmap_x(_day)
                    if _x_data is None:
                        continue
                    _x_disp = _clusterMap.ax_heatmap.transData.transform((_x_data, 0))[0]
                    _x_fig = _fig.transFigure.inverted().transform((_x_disp, 0))[0]
                    _fig.add_artist(_mlines.Line2D([_x_fig, _x_fig], [_y_bot, _y_top], transform=_fig.transFigure, color='black', linewidth=1.2, linestyle='--', alpha=0.6))
        if mode == 'log2fc':
            # Align the dendrogram's left edge with the right edge of the longest y-tick label.
            _fig = _clusterMap.figure
            _fig.canvas.draw()
            _renderer = _fig.canvas.get_renderer()
            _max_right_disp = None
            for _lab in _clusterMap.ax_heatmap.get_yticklabels():
                if not _lab.get_text():
                    continue
                _bb = _lab.get_window_extent(renderer=_renderer)
                if _max_right_disp is None or _bb.x1 > _max_right_disp:
                    _max_right_disp = _bb.x1
            if _max_right_disp is not None:
                _max_right_fig = _fig.transFigure.inverted().transform((_max_right_disp, 0))[0]
                _dend_pos = _clusterMap.ax_row_dendrogram.get_position()
                _clusterMap.ax_row_dendrogram.set_position([_max_right_fig + 0.008, _dend_pos.y0, dendro_w, _dend_pos.height])
        _suffix = '_diff' if mode == 'log2fc' else ''
        _clusterMap.figure.savefig(f"{title.lower().replace(' ', '_')}{_suffix}.png", bbox_inches='tight', dpi=300)
    from numpy import log, log2
    _DEFAULT_COLOR = 'lightgray'
    iterativeIDs_1 = json.load(open('iterativeIDs.json', 'r'))
    _zero_level = 1e-05

    # Shannon-diversity index computation disabled per request.
    # def shannon_index(abundances):
    #     return -sum([abundance * log(abundance) for abundance in abundances if abundance > _zero_level])
    # shannon_indices = {_sample: shannon_index(list(abs.to_numpy())) for _sample, abs in abundances.iterrows()}
    shannon_indices = None
    taxonomy_1 = json.load(open('iterativeID_taxonomy.json', 'r'))
    dic, _taxonomies = ({}, {})
    level_1 = 'Genus'
    for _sample, abs in abundances.iterrows():
        for _org, _ab in abs.items():
            if _ab <= _zero_level:
                continue
            _taxonomies.setdefault(_org, '|'.join([str(taxonomy_1[_org][l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index('Genus')]))
            dic.setdefault(_sample, {})
            dic[_sample].setdefault(_org, 0)
            dic[_sample][_org] = dic[_sample][_org] + _ab
    _nonzero_per_day = {_day: dict(sorted({_k: _v for _k, _v in _org_dict.items() if _v > 0}.items(), key=lambda item: item[1], reverse=True)) for _day, _org_dict in dic.items()}
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
        _top_per_day[_day] = {_org: log10(_v) for _org, _v in _org_dict.items() if _org in top10_all_days}
    _taxonomies = {_org: _taxa for _org, _taxa in _taxonomies.items() if _org in top10_all_days}
    df_1 = DataFrame(_top_per_day)
    df_1 = df_1.astype(float).replace([inf, -inf], nan)
    _taxonomy_series = Series({idx: _taxonomies.get(idx, f'Unknown|{idx}') for idx in df_1.index})
    mo.output.append(df_1)
    _arr = df_1.values.astype(float)
    row_means = nanmean(_arr, axis=1, keepdims=True)
    arr_filled = where(isnan(_arr), row_means, _arr)
    row_linkage = linkage(arr_filled, method='average', metric='euclidean')
    create_heatmap(df_1, _taxonomy_series, f'Top {_topNum} ASVs (% abundance)', shannon_indices, row_linkage)

    _root_abundances = abundances.copy()
    _root_abundances.columns = [c.split('.')[0] for c in _root_abundances.columns]
    _root_abundances = _root_abundances.T.groupby(level=0).sum().T
    _min_max_pct = 0.01
    _root_abundances = _root_abundances.loc[:, _root_abundances.max(axis=0) >= _min_max_pct]
    _root_taxonomies = {}
    for _ID, _taxa in taxonomy_1.items():
        _root = _ID.split('.')[0]
        if _root in _root_abundances.columns and _root not in _root_taxonomies:
            _root_taxonomies[_root] = '|'.join([str(_taxa[l]) for l in _taxonomic_levels if _taxonomic_levels.index(l) <= _taxonomic_levels.index('Genus')])

    # Per-phase mean relative abundance, then log2 fold change for each
    # unordered phase pair (all C(n_phases, 2) permutations explored).
    _phase_day_ranges = json.load(open('phase_day_ranges.json', 'r'))
    _phase_means = {}
    for _phase, _span in _phase_day_ranges.items():
        _start, _end = (_span['start'], _span['end'])
        _in_phase = [_d for _d in _root_abundances.index if _start <= int(_d) <= _end]
        if not _in_phase:
            continue
        _phase_means[_phase] = _root_abundances.loc[_in_phase].mean(axis=0)
    _phase_mean_df = DataFrame(_phase_means)
    _pseudocount = 1e-06
    _phase_order = [_p for _p in _phase_day_ranges.keys() if _p in _phase_mean_df.columns]
    _fc_cols = {}
    for _i, _p1 in enumerate(_phase_order):
        for _p2 in _phase_order[_i + 1:]:
            _label = f'{_p1} / {_p2}'
            _fc_cols[_label] = log2((_phase_mean_df[_p1] + _pseudocount) / (_phase_mean_df[_p2] + _pseudocount))
    df_root = DataFrame(_fc_cols)
    df_root = df_root.astype(float).replace([inf, -inf], nan)
    _root_taxonomy_series = Series({idx: _root_taxonomies.get(idx, f'Unknown|{idx}') for idx in df_root.index})
    _arr_root = df_root.values.astype(float)
    _row_means_root = nanmean(_arr_root, axis=1, keepdims=True)
    _arr_filled_root = where(isnan(_arr_root), _row_means_root, _arr_root)
    _row_linkage_root = linkage(_arr_filled_root, method='average', metric='euclidean')
    create_heatmap(df_root, _root_taxonomy_series, 'Taxa Above 1% Max Abundance', None, _row_linkage_root, mode='log2fc')
    return (
        Patch,
        Series,
        abs,
        colorsys,
        inf,
        iterativeID_color_map_1,
        json,
        nan,
        phylum_color_map,
        sns,
    )


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


@app.cell(disabled=True)
def _(DataFrame, Series, dump, load, spearmanr):
    from pandas import isna
    iterativeIDs_2 = load(open('model_inputs/iterativeIDs.json', 'r'))
    abundances_1 = load(open('model_inputs/abundances.json', 'r'))
    asvSet_abundances = load(open(f'modeling_files/ASVset_abundances.json', 'r'))
    omitted_columns = {'all_relevant_samples': ['10AB', 'A34', 'B34', 'C3', 'D34', 'E34', 'F34', 'G12', 'G3', 'H34']}

    def _correlations(ser1, ser2):
        aligned1, aligned2 = ser1.align(ser2, join='inner')
        if len(set(aligned1)) == 1 or len(set(aligned2)) == 1:
            global constant_vals
            constant_vals = constant_vals + 1
            return (float('nan'), float('nan'))
        return spearmanr(aligned1, aligned2)
    abundances_df = DataFrame(abundances_1).drop(ommitted, axis=1)
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
    return (
        ASVset_correlations,
        IterativeID_correlations,
        abundances_1,
        ommitted,
    )


@app.cell(disabled=True)
def _(ASVset_correlations, IterativeID_correlations, mo):
    mo.output.append(IterativeID_correlations)
    print(len(IterativeID_correlations), len(ASVset_correlations))
    return


@app.cell(disabled=True)
def _(DataFrame, Series, abs, load, mo, multipletests, read_csv, sns):
    from matplotlib import colors, patches
    from pandas import set_option
    from glob import glob
    set_option('display.max_rows', None)
    _taxonomic_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    _sample_days = {'10AB': 0, 'I34': 149, 'J12': 161, 'K12': 182, 'L12': 189, 'M12': 210, 'N12': 238, 'P12': 283, 'Q': 300}
    order = ['all samples $H_2$ feed', 'all samples $H_2$ BT']
    _correlations, _pvals = ({}, {})
    for _k in order:
        _correlations.setdefault(_k, {})
        _pvals.setdefault(_k, {})
    for cor in glob('modeling_files/correlations/ASVset_correlations_all_relevant_samples_H2*.json'):
        _name = cor.split('/')[-1].split('.')[0].split('_correlations_')[1].replace('_', ' ').replace('relevant ', '').replace('phase ', '').replace('H2', '$H_2$')
        content = load(open(cor, 'r'))
        _correlations[_name].update({_k: _v['correlation'] for _k, _v in content.items()})
        _pvals[_name].update({_k: _v['p_value'] for _k, _v in content.items()})
    q_vals = {}
    for _name, content in _pvals.items():
        pvals_list = list(content.values())
        q_vals[_name] = multipletests(pvals_list, alpha=0.05, method='fdr_bh')
        print(sum(q_vals[_name][0]))
    new_correlations = {_name: dict(inner) for _name, inner in _correlations.items()}
    new_pvals = {_name: dict(inner) for _name, inner in _pvals.items()}
    for _name, content in _pvals.items():
        for _i, (_k, _v) in enumerate(content.items()):
            if not q_vals[_name][0][_i]:
                new_correlations[_name].pop(_k)
                new_pvals[_name].pop(_k)
    mo.output.append(new_correlations)
    mo.output.append(new_pvals)
    df_2 = DataFrame(new_correlations).fillna(0)
    _pval_matrix = DataFrame(new_pvals)
    df_2.rename(columns={'all samples $H_2$ feed': '$H_2$ loading $\\left(\\frac{mol}{min}\\right)$', 'all samples $H_2$ BT': '$H_2$ breakthrough $\\left(\\frac{mol}{min}\\right)$'}, inplace=True)
    _pval_matrix.rename(columns={'all samples $H_2$ feed': '$H_2$ loading', 'all samples $H_2$ BT': '$H_2$ breakthrough'}, inplace=True)
    reduced = True
    if reduced:
        df_2 = df_2[(df_2 != float(0)).any(axis=1)]
    _total_df = read_csv('model_inputs/total.csv').set_index('seq')
    orgs = {}
    for _i in df_2.index:
        orgs.setdefault(_i.split('.')[0], []).append(_i)
    _taxonomies = {}
    level_2 = 'Genus'
    for _seq, _row in _total_df.iterrows():
        _day = _sample_days.get(_row['sample'])
        taxonomy_2 = []
        for l in reversed(_taxonomic_levels):
            _taxa = str(_row[l])
            IDs = orgs.get(_taxa)
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
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {_k.split('.')[0]: _v for _k, _v in _iterativeID_levels.items()}
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
    return df_2, orgs, patches


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


@app.cell(disabled=True)
def _(
    abundances_1,
    iterativeID_color_map_1,
    iterativeID_phylums_1,
    iterativeID_taxonomy_1,
):
    missing = [_i for _i in abundances_1.index if _i not in iterativeID_color_map_1]
    print('missing from color map:', missing)
    for _i in missing[:5]:
        print(_i, '→ taxonomy:', iterativeID_taxonomy_1.get(_i), 'phylum:', iterativeID_phylums_1.get(_i))
    return


@app.cell(disabled=True)
def _(
    DataFrame,
    GAOs_PAOs,
    Patch,
    Series,
    df_2,
    json,
    load,
    mo,
    nan,
    orgs,
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
    iterativeID_color_map_2 = json.load(open(f'iterativeID_color_map.json', 'r'))
    iterativeIDs_3 = json.load(open('iterativeIDs.json', 'r'))
    abundances_2 = read_csv('abundances.csv', header=0).set_index('sample')
    mo.output.append(abundances_2)
    abundances_2.drop([col for col in abundances_2.columns if col not in _significantly_connected_organisms], axis=1, inplace=True)
    abundances_2.drop([col for col in abundances_2.columns if abundances_2[col].max() < 0.005], axis=1, inplace=True)
    _total_captured = abundances_2.sum(axis=1)
    print(_total_captured)
    mo.output.append(abundances_2.head())
    _pval_matrix, corr_matrix = corr_pvalues(abundances_2)
    mo.output.append(corr_matrix.head())
    mo.output.append(_pval_matrix.head())
    taxonomy_3 = json.load(open('iterativeID_taxonomy.json', 'r'))
    _taxonomies = {}
    for col in abundances_2.columns:
        _taxonomies[col] = '|'.join([_v for _k, _v in taxonomy_3.get(col, 'Unknown').items() if _k != 'Species' and _v is not None])
    _taxonomy_series = Series({idx: _taxonomies.get(idx, f'Unknown|{idx}') for idx in abundances_2.columns})
    print(f"Detected depth: {max((len(t.split('|')) for t in _taxonomy_series))}")
    print(f'df rows: {len(df_2.columns)}')
    print(f'taxonomy_series length: {len(_taxonomy_series)}')
    print(f'Sample entries:\n{_taxonomy_series.head()}')
    taxa_color_map_1 = json.load(open('iterativeID_color_map.json'))
    _DEFAULT_COLOR = 'lightgray'
    bar_label = 'Phylum'
    row_colors = Series({idx: taxa_color_map_1.get(idx, _DEFAULT_COLOR) for idx in corr_matrix.index}, name=bar_label)
    col_colors = Series({col: taxa_color_map_1.get(col, _DEFAULT_COLOR) for col in corr_matrix.columns}, name=bar_label)
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
    _iterativeID_levels = json.load(open('iterativeID_levels.json', 'r'))
    _ID_levels = {}
    for _k, _v in _iterativeID_levels.items():
        _ID_levels[_k] = _v
        _ID_levels.setdefault(_k.split('.')[0], _v)
    _iterativeID_levels = json.load(open('iterativeID_levels.json', 'r'))
    inverted_GAOs_PAOs = {_v: _k for _k, vs in GAOs_PAOs.items() for _v in vs}
    for _label in _clusterMap.ax_heatmap.get_yticklabels():
        _text = _label.get_text()
        for _org, val in inverted_GAOs_PAOs.items():
            if _org not in _text:
                continue
            _label.set_fontweight('bold')
            if 'GAOs' in val:
                _label.set_color('green')
                if 'Putative' in val:
                    _label.set_color('mediumseagreen')
            elif 'PAOs' in val:
                _label.set_color('blue')
                if 'Putative' in val:
                    _label.set_color('cornflowerblue')
            elif 'PHA' in val:
                _label.set_color('red')
        if _ID_levels.get(_text) == 'Genus':
            _label.set_fontstyle('italic')
    for _label in _clusterMap.ax_heatmap.get_xticklabels():
        _text = _label.get_text()
        for _org, val in inverted_GAOs_PAOs.items():
            if _org not in _text:
                continue
            _label.set_fontweight('bold')
            if 'GAOs' in val:
                _label.set_color('green')
                if 'Putative' in val:
                    _label.set_color('mediumseagreen')
            elif 'PAOs' in val:
                _label.set_color('blue')
                if 'Putative' in val:
                    _label.set_color('cornflowerblue')
            elif 'PHA' in val:
                _label.set_color('red')
        if _ID_levels.get(_text) == 'Genus':
            _label.set_fontstyle('italic')
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
    for _org in orgs:
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
    _gao_pao_box_color = {'GAOs': 'green', 'Putative GAOs': 'mediumseagreen', 'PAOs': 'blue', 'Putative PAOs': 'cornflowerblue'}
    _n_corr = len(corr_matrix.index)
    for _id in corr_matrix.index:
        _cat = inverted_GAOs_PAOs.get(_id.split('.')[0])
        if _cat not in _gao_pao_box_color:
            continue
        _box_color = _gao_pao_box_color[_cat]
        _row_pos = _dendrogram_row.index(corr_matrix.index.get_loc(_id))
        _col_pos = _dendrogram_col.index(corr_matrix.columns.get_loc(_id))
        if _one_triangle:
            _row_rect = patches.Rectangle((0, _row_pos), _row_pos + 1, 1, linewidth=4, edgecolor=_box_color, facecolor='none', clip_on=False)
            _col_rect = patches.Rectangle((_col_pos, _n_corr), 1, -(_n_corr - _col_pos), linewidth=4, edgecolor=_box_color, facecolor='none', clip_on=False)
        else:
            _row_rect = patches.Rectangle((0, _row_pos), len(corr_matrix.columns), 1, linewidth=4, edgecolor=_box_color, facecolor='none', clip_on=False)
            _col_rect = patches.Rectangle((_col_pos, 0), 1, _n_corr, linewidth=4, edgecolor=_box_color, facecolor='none', clip_on=False)
        _clusterMap.ax_heatmap.add_patch(_row_rect)
        _clusterMap.ax_heatmap.add_patch(_col_rect)
    proteo_class_color_1 = json.load(open('proteo_class_color.json'))
    _proteo_base = json.load(open('phylum_base_overrides.json'))['Proteobacteria']
    if isinstance(row_colors, Series):
        rc = row_colors
    else:
        rc = Series(list(row_colors), index=df_2.index)
    phylum_color = {}
    for idx in corr_matrix.index:
        parts = str(_taxonomies.get(idx, '')).split('|')
        if len(parts) < 2:
            continue
        phylum = parts[1]
        if phylum in ('None', '', 'Unknown'):
            continue
        if phylum == 'Proteobacteria':
            phylum_color['Proteobacteria'] = _proteo_base
            continue
        _color = rc.get(idx)
        if _color is None:
            continue
        phylum_color.setdefault(phylum, _color)
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
            if _p == 'Proteobacteria' and proteo_class_color_1:
                handles.append(Patch(color='none', label=_p))
                for cls, _color in proteo_class_color_1.items():
                    handles.append(Patch(facecolor=_color, label=f'      {cls}'))
            else:
                handles.append(Patch(facecolor=phylum_color[_p], label=_p))
    _clusterMap.figure.legend(handles=handles, title='Phylum', title_fontsize=50, fontsize=40, loc='upper right', bbox_to_anchor=(0.65, 0.7), frameon=True, borderaxespad=0.5, handlelength=1.5, handletextpad=0.6)
    _clusterMap.ax_heatmap.set_xlabel('Member ASVs', fontsize=50)
    _clusterMap.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)
    _clusterMap.figure.savefig(f"abundance_correlatons{('_one_triangle' if _one_triangle else '')}.png", dpi=300, bbox_inches='tight')
    return corr_matrix, ones, ones_like, proteo_class_color_1, triu


@app.cell(disabled=True)
def _(proteo_class_color_1):
    proteo_class_color_1
    return


@app.cell(disabled=True)
def _(KMeans, best_k, corr_matrix, pd, sig_mask):
    # sig_mask is a boolean DataFrame (True = p < threshold)
    corr_masked = corr_matrix.where(sig_mask, 0)
    km = KMeans(n_clusters=best_k, n_init=20, random_state=42)
    cluster_labels = pd.Series(km.fit_predict(corr_masked), index=corr_matrix.index)
    cluster_labels
    return


@app.cell(disabled=True)
def _(
    DataFrame,
    abs,
    inf,
    load,
    mo,
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
        orgs_1 = dict(list(_org_dict.items())[:_topNum])
        _top_per_day[_day] = orgs_1
        _all_orgs.update(orgs_1)
    df_3 = DataFrame(_top_per_day).T
    df_3 = df_3.astype(float).replace([inf, -inf], nan)
    df_3 = df_3.loc[:, df_3.notna().sum() >= 5]
    mo.output.append(df_3)
    corr_matrix_1 = df_3.corr('spearman').dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    mo.output.append(corr_matrix_1)
    _pval_matrix = corr_pvalues_1(df_3).dropna(axis=1, how='all').dropna(axis=0, how='all').fillna(0)
    _pval_matrix = _pval_matrix[corr_matrix_1.columns]
    _pval_matrix = _pval_matrix.loc[corr_matrix_1.index]
    mo.output.append(_pval_matrix)
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
    _iterativeID_levels = load(open('model_inputs/iterativeID_levels.json', 'r'))
    _ID_levels = {_k.split('.')[0]: _v for _k, _v in _iterativeID_levels.items()}
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
    xlabels = _clusterMap.ax_heatmap.get_xticklabels()
    for _label in xlabels:
        if _label.get_text() in _organisms_to_highlight:
            _label.set_fontsize(_labelsize * 1.2)
            _label.set_fontweight('bold')
        _text = _label.get_text().split('.')[0]
        _taxa = _ID_levels.get(_text)
        if _taxa == 'Genus':
            _label.set_fontstyle('italic')
    _clusterMap.ax_heatmap.set_xticklabels(xlabels)
    _clusterMap.ax_row_dendrogram.set_visible(False)
    _clusterMap.ax_col_dendrogram.set_visible(False)
    _heatmap_pos = _clusterMap.ax_heatmap.get_position()
    _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 1.01, _heatmap_pos.y0 + 0.5, 0.06, _heatmap_pos.height / 5])
    if _one_triangle:
        _clusterMap.ax_cbar.set_position([_heatmap_pos.x1 - 0.3, _heatmap_pos.y0 + 0.4, 0.06, _heatmap_pos.height / 5])
    _clusterMap.ax_heatmap.set_xlabel('ASVs', fontsize=100, labelpad=20)
    _clusterMap.ax_heatmap.set_ylabel('ASVs', fontsize=100, labelpad=20)
    _clusterMap.figure.savefig(f"abundance_heatmaps/Top_{_topNum}_ASVs_abundance_correlation{('_one_triangle' if _one_triangle else '')}.png", bbox_inches='tight', dpi=300)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # co-occurrence figure
    """)
    return


@app.cell
def _(
    DataFrame,
    GAOs_PAOs,
    combinations,
    defaultdict,
    json,
    load,
    mcolors,
    mpatches,
    np,
    pair_data,
    plt,
    read_csv,
    reject,
):
    import networkx as nx
    import matplotlib.cm as cm
    from numpy import ones as _np_ones
    from scipy.stats import spearmanr as _spearmanr_local

    iterativeID_taxonomy_local = load(open('iterativeID_taxonomy.json'))
    level_3 = 'Phylum'
    iterativeID_level = {_ID: content.get(level_3, 'Unknown') for _ID, content in iterativeID_taxonomy_local.items()}

    _significantly_connected_organisms = [str(x) for x in np.load('FDR_passing_pairs.npy')]
    _significantly_connected_organisms.append('Methanobacteriaceae.1')

    iterativeID_color_map_local = json.load(open('iterativeID_color_map.json', 'r'))
    taxa_color_map_local = iterativeID_color_map_local
    order_color_map_local = json.load(open('Phylum_color_map.json', 'r'))

    abund_thresh = 0.005  # 0.5% mean relative abundance
    abundances_3 = read_csv('abundances.csv', header=0).set_index('sample')
    abundances_3.drop([col for col in abundances_3.columns if col not in _significantly_connected_organisms], axis=1, inplace=True)

    def _corr_pvalues_local(df):
        n = df.shape[1]
        pv = DataFrame(_np_ones((n, n)), index=df.columns, columns=df.columns)
        cr = DataFrame(_np_ones((n, n)), index=df.columns, columns=df.columns)
        for _i in range(n):
            for _j in range(_i + 1, n):
                c, _p = _spearmanr_local(df.iloc[:, _i], df.iloc[:, _j])
                pv.iloc[_i, _j] = _p; pv.iloc[_j, _i] = _p
                cr.iloc[_i, _j] = c; cr.iloc[_j, _i] = c
        return (pv, cr)

    pvals_matrix, corr_matrix_2 = _corr_pvalues_local(abundances_3)
    _mean_rel_abund = abundances_3.mean(axis=0)

    _presence = (abundances_3 > 0).astype(int)
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
        import math as _math
        import matplotlib.patheffects as path_effects
        from scipy.spatial import ConvexHull as _ConvexHull

        # === Compute Louvain communities (positive-rho subgraph) + singleton handling for isolates ===
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
            print(f'[{title_slug}] isolates: {len(_isolate_nodes)} nodes with only negative-rho edges -> singleton modules')

        _multi_comms_only = [c for c in _layout_comms if len(c) >= 2]
        _n_singletons = sum(1 for c in _layout_comms if len(c) == 1)

        # === Module-membership JSON export ===
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
        print(f'[{title_slug}] wrote network_module_membership_{title_slug}.json ({_louvain_n} louvain modules + {_singleton_n} singletons = {len(_module_export)} total)')

        # === Anchor warm-start: place each multi-node module at a circular anchor ===
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

        # === Layout: per-module subgraph layout + grid/ring placement, or unweighted spring_layout ===
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
                # === Grid placement: rotated sqrt(n) x sqrt(n) grid with vertical stretch ===
                _intra_amp = 5.4 * 3
                for _mi in _module_layouts:
                    for _n in list(_module_layouts[_mi].keys()):
                        _module_layouts[_mi][_n] = _module_layouts[_mi][_n] * _intra_amp
                _max_module_radius = max(
                    (max(np.linalg.norm(p) for p in _module_layouts[_mi].values())
                     for _mi in _module_layouts), default=1.0)
                _hull_buffer = 1.0
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
                    _singleton_radius = max(_xext, _yext) + _max_module_radius * 0.2
                    for _i, _mi in enumerate(_singleton_mis):
                        _theta = 2 * np.pi * (_i + 0.5) / max(len(_singleton_mis), 1)
                        for _n in _layout_comms[_mi]:
                            pos[_n] = np.array([np.cos(_theta), np.sin(_theta)]) * _singleton_radius
            else:
                # === Ring placement: modules around a circle (default) ===
                _max_module_radius = max(
                    (max(np.linalg.norm(p) for p in _module_layouts[_mi].values())
                     for _mi in _module_layouts), default=1.0)
                _hull_buffer = 1.02
                _gap_factor = 1.01
                _n_multi_mods = len(_module_layouts)
                _min_chord = 2 * _max_module_radius * _hull_buffer * _gap_factor
                _sep_radius = _min_chord / (2 * np.sin(np.pi / _n_multi_mods)) if _n_multi_mods > 1 else 0
                print(f'[{title_slug}] per-module layout: ring radius={_sep_radius:.2f}, max module radius={_max_module_radius:.2f}')
                pos = {}
                _sorted_mis = sorted(_module_layouts.keys(), key=lambda mi: -len(_layout_comms[mi]))
                for _i, _mi in enumerate(_sorted_mis):
                    _theta = 2 * np.pi * _i / max(_n_multi_mods, 1)
                    _module_center = np.array([np.cos(_theta), np.sin(_theta)]) * _sep_radius
                    for _n, _local in _module_layouts[_mi].items():
                        pos[_n] = _module_center + _local
                _singleton_mis = [_mi for _mi, _c in enumerate(_layout_comms) if len(_c) == 1]
                if _singleton_mis:
                    _singleton_radius = _sep_radius + _max_module_radius * 2.5
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
            _module_centroid = {}
            for _mi, _comm in enumerate(_layout_comms):
                _pts = np.array([pos[_n] for _n in _comm if _n in pos])
                if len(_pts) == 0:
                    continue
                _module_centroid[_mi] = _pts.mean(axis=0)
            _pair_pos_sum = defaultdict(float)
            _pair_pos_count = defaultdict(int)
            _pair_neg_sum = defaultdict(float)
            _pair_neg_count = defaultdict(int)
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
            _count_legend_dy = 0.035
            _count_legend_top_y = 4 * _count_legend_dy + 0.01  # 0.15 — bottom of thickest line at ~y=0
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

        # --- Module shading: SAT-shrunk convex hulls ---
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
        _initial_scale = 1.30
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

        _min_scale = 0.75
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
        plt.savefig(f'cooccurrence_network_{title_slug}_min{abund_thresh*100:g}pct{_layout_suffix}{_gamma_suffix}{_simplified_suffix}.png', dpi=300, bbox_inches='tight')
        plt.show()
        plt.close(fig)

    genus_to_ids = defaultdict(list)
    for _ID, _taxa in iterativeID_taxonomy_local.items():
        g = _taxa.get('Genus', '')
        if g:
            genus_to_ids[g].append(_ID)

    def resolve_to_graph_nodes(genus_names, graph):
        ids = []
        for genus in genus_names:
            for _ID in genus_to_ids.get(genus, []):
                if _ID in graph:
                    ids.append(_ID)
        return ids

    # Full-network grid-layout render (emulates codiffusion_bioreactor's cooccurrence_network_p_value_FDR_grid_layout.png)
    render_network(G, 'p_value_FDR', grid_layout_mode=True, resolution=1.115)

    selections = {'GAOs': GAOs_PAOs['GAOs'], 'PAOs': GAOs_PAOs['PAOs']}
    for _name, genus_names in selections.items():
        focus_in_G = resolve_to_graph_nodes(genus_names, G)
        print(f'[{_name}] genera requested: {genus_names}')
        print(f'[{_name}] resolved {len(focus_in_G)} iterativeIDs in G: {focus_in_G}')
        if not focus_in_G:
            print(f'[{_name}] no focus nodes are in G — skipping')
            continue
        focus_set = set(focus_in_G)
        high_abund = {_n for _n in G.nodes() if _mean_rel_abund.get(_n, 0) > abund_thresh}
        keepable = focus_set | high_abund
        incident_edges = [(u, v) for u, v in G.edges()
                          if (u in focus_set or v in focus_set) and u in keepable and v in keepable]
        G_sub = G.edge_subgraph(incident_edges).copy()
        render_network(G_sub, _name, focus_nodes=set(focus_in_G))

    return


@app.cell(disabled=True)
def _(GAOs_PAOs):
    """Per-phase versions of the one-triangle correlation and co-occurrence network figures."""
    from json import load as _load
    from collections import defaultdict as _dd
    from itertools import combinations as _comb
    import numpy as _np
    from numpy import triu as _triu, ones_like as _ones_like, nan as _nan, ones as _ones
    from pandas import DataFrame as _DF, Series as _Ser, read_csv as _read_csv
    from scipy.stats import spearmanr as _spr
    from statsmodels.stats.multitest import multipletests as _mt
    import matplotlib.pyplot as _plt
    import matplotlib.colors as _mc
    import matplotlib.cm as _cm_mpl
    import matplotlib.patches as _mpatch
    import seaborn as _sns
    import networkx as _nx
    _sd = _load(open('sample_days.json'))
    _phr = _load(open('phase_day_ranges.json'))
    _itx = _load(open('iterativeID_taxonomy.json'))
    _icm = _load(open('iterativeID_color_map.json'))
    _ilv = _load(open('iterativeID_levels.json'))
    _pcm = _load(open('Phylum_color_map.json'))
    _pcc = _load(open('proteo_class_color.json'))
    _pb = _load(open('phylum_base_overrides.json'))['Proteobacteria']
    _inv = {v: k for k, vs in GAOs_PAOs.items() for v in vs}
    _LBL = {'GAOs': 'green', 'Putative GAOs': 'mediumseagreen', 'PAOs': 'blue', 'Putative PAOs': 'cornflowerblue', 'Other PHA storing potential+ function': 'red'}
    _BOX = {'GAOs': 'green', 'Putative GAOs': 'mediumseagreen', 'PAOs': 'blue', 'Putative PAOs': 'cornflowerblue'}
    _IDLV = {}
    for _k, _v in _ilv.items():
        _IDLV[_k] = _v
        _IDLV.setdefault(_k.split('.')[0], _v)
    _ABT = 0.005
    _ZL = 1e-05
    _EXCLUDE_DAYS_BELOW = 15

    def _gpc(node_id, mapping):
        text = str(node_id)
        for _o, _c in _inv.items():
            if _o in text:
                return mapping.get(_c)
        return None

    def _corr_pv(df):
        n = df.shape[1]
        pv = _DF(_ones((n, n)), index=df.columns, columns=df.columns)
        cr = _DF(_ones((n, n)), index=df.columns, columns=df.columns)
        for _i in range(n):
            for _j in range(_i + 1, n):
                _c, _p = _spr(df.iloc[:, _i], df.iloc[:, _j])
                pv.iloc[_i, _j] = _p; pv.iloc[_j, _i] = _p
                cr.iloc[_i, _j] = _c; cr.iloc[_j, _i] = _c
        return (pv, cr)

    def _fdr(df):
        pres = (df > _ZL).astype(int)
        cooc = _dd(int)
        for s in pres.itertuples(index=False):
            present = [c for c, v in zip(pres.columns, s) if v]
            for p in _comb(sorted(present), 2):
                cooc[p] += 1
        pdat = []
        for (_a, _b), _ct in cooc.items():
            _r, _pv = _spr(df[_a], df[_b])
            if _np.isnan(_r):
                continue
            pdat.append((_a, _b, _r, _pv, _ct))
        if not pdat:
            return [], set()
        pvs = _np.array([t[3] for t in pdat])
        rj, _, _, _ = _mt(pvs, alpha=0.05, method='fdr_bh')
        passing = [pdat[_i] for _i in range(len(pdat)) if rj[_i]]
        orgs = {a for a, b, *_ in passing} | {b for a, b, *_ in passing}
        return passing, orgs

    def _build_tax(orgs):
        return _Ser({o: '|'.join([_itx.get(o, {}).get(l, '') for l in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']]) for o in orgs})

    def _render_one_triangle_phase(corr_matrix, taxonomies, suffix):
        if corr_matrix.shape[0] < 2:
            print(f'[one-triangle {suffix}] not enough orgs'); return
        DEFAULT = 'lightgray'
        rc = _Ser({i: _icm.get(i, DEFAULT) for i in corr_matrix.index}, name='Phylum')
        cc = _Ser({c: _icm.get(c, DEFAULT) for c in corr_matrix.columns}, name='Phylum')
        cm = _sns.clustermap(corr_matrix, row_colors=rc, col_colors=cc, cmap='coolwarm_r', center=0, figsize=(60, 70), dendrogram_ratio=(0.1, 0.2))
        cm.figure.subplots_adjust(bottom=0.15, top=0.95)
        cm.ax_row_dendrogram.set_visible(False); cm.ax_col_dendrogram.set_visible(False)
        cm.ax_heatmap.yaxis.set_ticks_position('left'); cm.ax_heatmap.yaxis.set_label_position('left')
        d_row = cm.dendrogram_row.reordered_ind; d_col = cm.dendrogram_col.reordered_ind
        df_re = corr_matrix.iloc[d_row, d_col]
        mask = _triu(_ones_like(df_re, dtype=bool), k=1)
        mesh = cm.ax_heatmap.collections[0]
        arr = mesh.get_array().reshape(df_re.shape); arr[mask] = _nan; mesh.set_array(arr.ravel())
        for axis_get in (cm.ax_heatmap.get_yticklabels, cm.ax_heatmap.get_xticklabels):
            for lbl in axis_get():
                t = lbl.get_text()
                for _o, _c in _inv.items():
                    if _o in t:
                        lbl.set_fontweight('bold')
                        col = _LBL.get(_c)
                        if col:
                            lbl.set_color(col)
                if _IDLV.get(t) == 'Genus':
                    lbl.set_fontstyle('italic')
        N = len(corr_matrix.index)
        for _id in corr_matrix.index:
            cat = _inv.get(_id.split('.')[0])
            if cat not in _BOX:
                continue
            color = _BOX[cat]
            rp = d_row.index(corr_matrix.index.get_loc(_id))
            cp = d_col.index(corr_matrix.columns.get_loc(_id))
            cm.ax_heatmap.add_patch(_mpatch.Rectangle((0, rp), rp + 1, 1, linewidth=4, edgecolor=color, facecolor='none', clip_on=False))
            cm.ax_heatmap.add_patch(_mpatch.Rectangle((cp, N), 1, -(N - cp), linewidth=4, edgecolor=color, facecolor='none', clip_on=False))
        cm.ax_heatmap.set_xlabel('Member ASVs', fontsize=50); cm.ax_heatmap.set_ylabel('Member ASVs', fontsize=50)
        out = f'abundance_correlatons_one_triangle_{suffix}.png'
        cm.figure.savefig(out, dpi=200, bbox_inches='tight'); _plt.close(cm.figure)
        print(f'wrote {out}')

    def _render_network_phase(G_sub, mean_rel, suffix, focus=None):
        if G_sub.number_of_edges() == 0:
            print(f'[network {suffix}] no edges'); return
        pos = _nx.spring_layout(G_sub, seed=42, iterations=100, k=0.3)
        edges = G_sub.edges(data=True)
        rho_values = [d['rho'] for _, _, d in edges]
        edge_widths = [3 * d['weight'] for _, _, d in edges]
        norm = _mc.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        cmap = _cm_mpl.RdBu
        edge_colors = [cmap(norm(r)) for r in rho_values]
        width = 40
        fig, ax = _plt.subplots(figsize=(width, 30))
        scale = 5000 * (width / 10) * 2
        cps = _np.sqrt
        node_sizes = [scale * cps(mean_rel.get(n, 0)) for n in G_sub.nodes()]
        node_colors = [_pcm.get(_itx.get(n, {}).get('Phylum', ''), 'lightgray') for n in G_sub.nodes()]
        edgecolors = [(_gpc(n, _LBL) or 'none') if focus and n in focus else 'none' for n in G_sub.nodes()]
        linewidths = [6 if focus and n in focus else 0 for n in G_sub.nodes()]
        _nx.draw_networkx_nodes(G_sub, pos, node_size=node_sizes, node_color=node_colors, edgecolors=edgecolors, linewidths=linewidths, alpha=0.9, ax=ax)
        _nx.draw_networkx_edges(G_sub, pos, width=edge_widths, edge_color=edge_colors, alpha=0.85, ax=ax)
        max_a = mean_rel.max() if len(mean_rel) else 0
        for n, (x, y) in pos.items():
            ab = mean_rel.get(n, 0)
            fs = 6 + (14 * cps(ab) / cps(max_a) * (width / 10) if max_a > 0 else 0)
            fs = max(5, min(fs, 48))
            color = _gpc(n, _LBL) or 'black'
            ax.text(x, y, str(n), fontsize=fs, color=color, fontweight='bold', ha='center', va='center')
        sm = _cm_mpl.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        cb = _plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cb.set_label('Spearman ρ', fontsize=10 * (width / 10))
        ax.axis('off'); _plt.tight_layout()
        out = f'cooccurrence_network_{suffix}.png'
        _plt.savefig(out, dpi=300, bbox_inches='tight'); _plt.close(fig)
        print(f'wrote {out}')

    _g2i = _dd(list)
    for _ID, _t in _itx.items():
        _g = _t.get('Genus', '')
        if _g:
            _g2i[_g].append(_ID)

    _abf = _read_csv('abundances.csv', header=0).set_index('sample')
    for _phase, _span in _phr.items():
        _lo = max(_span['start'], _EXCLUDE_DAYS_BELOW)
        print(f"\n=== Phase {_phase} (days {_lo}-{_span['end']}; first {_EXCLUDE_DAYS_BELOW - 1} days excluded) ===")
        _kp = [s for s in _abf.index if _lo <= int(_sd[s]) <= _span['end']]
        _ab = _abf.loc[_kp]
        _ab = _ab.loc[:, (_ab.fillna(0) > 0).any(axis=0)]
        print(f"  samples: {_ab.shape[0]}, organisms: {_ab.shape[1]}")
        if _ab.shape[0] < 3 or _ab.shape[1] < 2:
            continue
        _pp, _po = _fdr(_ab)
        if not _po:
            continue
        _absig = _ab[[c for c in _ab.columns if c in _po]]
        _abc = _absig.loc[:, _absig.max() >= _ABT]
        if _abc.shape[1] >= 2:
            _, _corr = _corr_pv(_abc)
            _render_one_triangle_phase(_corr, _build_tax(_corr.columns), f'phase{_phase}')
        _G = _nx.Graph()
        for _a, _b, _r, _pv, _ct in _pp:
            _G.add_edge(_a, _b, weight=abs(_r), rho=_r, pvalue=_pv, cooccurrence=_ct)
        _mr = _ab.mean(axis=0)
        for _sn, _gn in {'GAOs': GAOs_PAOs['GAOs'], 'PAOs': GAOs_PAOs['PAOs']}.items():
            _focus = [_id for _g in _gn for _id in _g2i.get(_g, []) if _id in _G]
            if not _focus:
                continue
            _high = {n for n in _G.nodes() if _mr.get(n, 0) > _ABT}
            _kpb = set(_focus) | _high
            _es = [(u, v) for u, v in _G.edges() if (u in set(_focus) or v in set(_focus)) and u in _kpb and v in _kpb]
            _Gs = _G.edge_subgraph(_es).copy()
            _render_network_phase(_Gs, _mr, f'{_sn}_phase{_phase}_min{_ABT*100:g}pct', focus=set(_focus))
@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Diff / innoculum / time / FDR heatmaps (isolated)

    This cell only reads pre-computed JSON/CSV artifacts and writes the four
    log2-fold-change heatmap PNGs. It has no marimo inputs, so re-running it
    will NOT trigger any of the upstream pipeline cells.
    """)
    return


@app.cell
def _():
    """Re-render the diff / innoculum / time / FDR heatmaps from disk.

    All figure logic lives in render_diff_heatmaps.py. The reload picks up
    any edits made there since the last run without restarting marimo.
    """
    import importlib
    import render_diff_heatmaps
    importlib.reload(render_diff_heatmaps)
    render_diff_heatmaps.main()
    return


if __name__ == "__main__":
    app.run()
