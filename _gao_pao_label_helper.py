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
