#!/usr/bin/env python3
"""
Bundle the irreducible set of files needed to iterate on the figure offline.

Contents of offline_bundle/:
  gene_rpkm.tsv            - RPKM matrix (KO groups × MAGs)
  gene_rpkm_per_sample.tsv - per-sample RPKM matrix
  selected_mags.tsv        - MAG list + classification + abundance
  taxonomy_labels.tsv      - Genus/species/iterativeID per MAG
  target_genes.bed         - target KO genes with genomic coords (reference only)
  04_plot_heatmap.py       - standalone plotting script
  README.md                - instructions
"""
import shutil, os

WORK = '/scratch1/afreiburger/emilykin/gene_ab_figure'
DATA = f'{WORK}/data'
BUNDLE = f'{WORK}/offline_bundle'

os.makedirs(BUNDLE, exist_ok=True)
os.makedirs(f'{BUNDLE}/data', exist_ok=True)

# Data files go into data/ subdir (matches what 04_plot_heatmap.py expects)
DATA_FILES = [
    (f'{DATA}/gene_rpkm.tsv',            'data/gene_rpkm.tsv'),
    (f'{DATA}/gene_rpkm_per_sample.tsv', 'data/gene_rpkm_per_sample.tsv'),
    (f'{DATA}/selected_mags.tsv',        'data/selected_mags.tsv'),
    (f'{DATA}/taxonomy_labels.tsv',      'data/taxonomy_labels.tsv'),
    (f'{DATA}/target_genes.bed',         'data/target_genes.bed'),
]
# Scripts go in root
SCRIPT_FILES = [
    (f'{WORK}/04_plot_heatmap.py',       '04_plot_heatmap.py'),
]
FILES = DATA_FILES + SCRIPT_FILES

for src, dst in FILES:
    if os.path.exists(src):
        shutil.copy2(src, f'{BUNDLE}/{dst}')
        print(f'  Copied: {dst}')
    else:
        print(f'  MISSING (skip): {src}')

# Write README
readme = """# Gene Abundance Figure — Offline Bundle

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
"""

with open(f'{BUNDLE}/README.md', 'w') as f:
    f.write(readme)
print('  Wrote: README.md')

# Patch the offline plot script to use the bundle directory
plot_src = f'{BUNDLE}/04_plot_heatmap.py'
with open(plot_src) as f:
    src = f.read()
# Replace the hardcoded WORK path with the local directory detection
patched = src.replace(
    "WORK = os.path.dirname(os.path.abspath(__file__))",
    "WORK = os.path.dirname(os.path.abspath(__file__))  # auto-set to this file's dir"
)
with open(plot_src, 'w') as f:
    f.write(patched)

# Create a zip (walk subdirectories)
import zipfile, pathlib
zip_path = f'{WORK}/offline_bundle.zip'
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
    for path in pathlib.Path(BUNDLE).rglob('*'):
        if path.is_file():
            zf.write(path, path.relative_to(BUNDLE))
print(f'\nCreated archive: {zip_path}')
print(f'  Size: {os.path.getsize(zip_path) / 1024:.0f} KB')

print('\nOffline bundle contents:')
for path in sorted(pathlib.Path(BUNDLE).rglob('*')):
    if path.is_file():
        rel = path.relative_to(BUNDLE)
        print(f'  {str(rel):40s}  {path.stat().st_size:>10,} bytes')
