# Gene Abundance Figure — Offline Bundle

## Contents
- `gene_rpkm.tsv`           — RPKM matrix: rows = KO groups, cols = MAG IDs
- `gene_rpkm_per_sample.tsv`— same but per CAN sample
- `selected_mags.tsv`       — selected MAGs: classification + relative abundance
- `taxonomy_labels.tsv`     — Genus/Species/iterativeID per MAG
- `target_genes.bed`        — target gene genomic coordinates (reference)
- `04_plot_heatmap.py`      — standalone plotting script

## How to regenerate the figure
```bash
# Edit 04_plot_heatmap.py to point WORK to this directory:
#   WORK = '/path/to/offline_bundle'
python3 04_plot_heatmap.py
```
Output: `gene_abundance_figure.png` and `gene_abundance_figure.pdf`

## Key customization points in 04_plot_heatmap.py

### Row labels
```python
ROWS = [
    ('nirB/D',    'label text', '#color'),
    ...
]
```

### Column ordering / filtering
Modify `sort_key()` or filter `sorted_mags` list.

### Classification colors
```python
CAT_COLORS = {'GAO': '#E6811A', 'PAO': '#2E86C1', ...}
```

### Value transformation
Change `np.log10(matrix + 1)` for linear scale:
```python
# Linear:
display_matrix = matrix
vmax = matrix.max()
# Log10:
display_matrix = np.log10(matrix + 1)
vmax = display_matrix.max()
```

### Add per-sample view
Load `gene_rpkm_per_sample.tsv` instead. Each MAG will have 5 sub-columns.

## Classification logic
MAGs are classified using flag columns from mag_abundance_summary.tsv:
- PAO: `stronger_putative_PAO_like = 1`
- GAO: `GAO_like = 1`
- Denitrifier: `complete_denitrification = 1` OR `incomplete_denitrification = 1`
Multiple flags = combined label (e.g., "PAO/Denitrifier").

## Dependencies
- Python ≥ 3.8
- matplotlib ≥ 3.5
- numpy ≥ 1.20
(No seaborn or scipy required)
