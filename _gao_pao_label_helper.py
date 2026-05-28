"""Shared GAO/PAO axis-label coloring for the operational correlation heatmaps.

Mirrors the substring-match-then-color logic used by the heatmap cells in
data_processing.py and the dbRDA species labels in dbRDA_analysis.py so every
figure uses one source of truth for category membership and colors.
"""
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
_INVERTED = {v: k for k, vs in GAOs_PAOs.items() for v in vs}

# ---------------------------------------------------------------------------
# Operational-variable selection for the correlation figures: the dbRDA
# variable set plus Ax_time. Split into performance (KPI outcomes, drawn above
# the dashed separator) and environmental/operational-driver rows (below).
# Parameter strings match the raw column names in the correlation tables.
PERFORMANCE_PARAMS = [
    'peakN2O [mg/L]',
    'specific denitrification rates [mg NO2–N g−1\xa0VSS−1 h−1]',
    'N removal (ppm) [N-ppn]',
    'P removal [P%]',
]
ENVIRONMENTAL_PARAMS = [
    'COD:N',
    'N:P',
    'DO_Avg_1350_1410',
    'DO_Max_1350_1410',
    'A_time [minutes]',
    'Ax_time [minutes]',
    'Acetate ppm [mg/L]',
    'Propionate ppm [mg/L]',
    'N_Ax-1 [mg/L]',
    'N_Ax-2 [mg/L]',
]
INCLUDED_PARAMS = PERFORMANCE_PARAMS + ENVIRONMENTAL_PARAMS

# Publication display labels (with units) for the heatmap parameter rows.
# Keyed by the raw column name in the correlation tables.
DISPLAY_LABELS = {
    # performance (top half)
    'specific denitrification rates [mg NO2–N g−1\xa0VSS−1 h−1]':
        'specific denitrification (mg NO₂⁻-N g⁻¹ VSS h⁻¹)',
    'N removal (ppm) [N-ppn]': 'N removal (NO₂⁻-N mg/L)',
    'P removal [P%]':          'P removal (P%)',
    'peakN2O [mg/L]':          'peak N₂O (mg N/L)',
    # environmental / operational drivers (bottom half)
    'DO_Avg_1350_1410':        'DO avg (mg O₂ L⁻¹)',
    'DO_Max_1350_1410':        'DO max (mg O₂ L⁻¹)',
    'A_time [minutes]':        'cumulative aeration (mg O₂ L⁻¹ cycle⁻¹)',
    'N_Ax-1 [mg/L]':           'N_Ax-1 (mg N/L)',
    'N_Ax-2 [mg/L]':           'N_Ax-2 (mg N/L)',
    'COD:N':                   'COD:N',
    'N:P':                     'N:P',
    'Acetate ppm [mg/L]':      'Acetate (mg/L)',
    'Propionate ppm [mg/L]':   'Propionate (mg/L)',
    'Ax_time [minutes]':       'Ax_time (min)',
}


def label_color_for(text: str) -> tuple[str, bool]:
    """Return (color, is_gao_pao) for a tick label text."""
    for org, cat in _INVERTED.items():
        if org in str(text):
            return GAO_PAO_LABEL_COLORS.get(cat, 'black'), True
    return 'black', False


def color_axis_labels(ax, axis: str = 'x') -> None:
    """Color and bold the GAO/PAO tick labels on the given axis ('x' or 'y')."""
    labels = ax.get_xticklabels() if axis == 'x' else ax.get_yticklabels()
    for label in labels:
        col, hit = label_color_for(label.get_text())
        if hit:
            label.set_color(col)
            label.set_fontweight('bold')


def cluster_order(mat, by_axis: int):
    """Leaf order for hierarchical clustering of rows (by_axis=0) or cols (by_axis=1)."""
    from scipy.cluster.hierarchy import linkage, leaves_list
    if mat.shape[by_axis] < 2:
        return list(range(mat.shape[by_axis]))
    X = mat.fillna(0).values
    if by_axis == 1:
        X = X.T
    Z = linkage(X, method='average', metric='euclidean')
    return list(leaves_list(Z))


def order_param_rows(rho_mat):
    """Filter a (parameter × organism) matrix to INCLUDED_PARAMS and order rows:
    performance KPI rows (clustered) on top, environmental/driver rows (clustered)
    below. Returns (reordered_matrix, n_performance_rows) so callers can place the
    dashed separator at row index = n_performance_rows.
    """
    rho_mat = rho_mat.loc[[p for p in rho_mat.index if p in INCLUDED_PARAMS]]
    perf_rows = [p for p in rho_mat.index if p in PERFORMANCE_PARAMS]
    env_rows = [p for p in rho_mat.index if p in ENVIRONMENTAL_PARAMS]
    if len(perf_rows) > 1:
        perf_rows = [perf_rows[i] for i in cluster_order(rho_mat.loc[perf_rows], 0)]
    if len(env_rows) > 1:
        env_rows = [env_rows[i] for i in cluster_order(rho_mat.loc[env_rows], 0)]
    out = rho_mat.reindex(perf_rows + env_rows)
    out = out.rename(index=DISPLAY_LABELS)  # raw column names → unit-bearing display labels
    return out, len(perf_rows)
