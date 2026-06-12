#!/usr/bin/env python3
"""
Gene abundance heatmap — denitrification KOs across top HQ MAGs.

Mirrors example_gene_ab.png layout:
  Rows  = KO functional groups
  Cols  = MAGs (union of top-10 HQ per sample)
  Cells = RPKM (log10-scaled display, raw RPKM annotated)
  Top bar = GAO / PAO / Denitrifier classification color

Requires: matplotlib, numpy, csv (stdlib)
No seaborn required — works offline.
"""

import csv, math, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.colorbar import ColorbarBase
import numpy as np

# ── Paths (adjust for offline use) ──────────────────────────────────────────
WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')

# ── Row definitions ──────────────────────────────────────────────────────────
# Colors match KO_groups.png annotation scheme:
#   green  = nitrite reductases (nirB/D DNRA + nirS denitrification)
#   salmon = NO / N₂O reductases (nirK + norC/nosZ)
ROWS = [
    ('nirB/D',    'nirB/D  –  NADH nitrite reductase\n(DNRA: NO₂ → NH₄)',   '#A5D6A7'),
    ('nirS',      'nirS  –  cytochrome cd₁ nitrite reductase\n(NO₂ → NO)',   '#A5D6A7'),
    ('nirK',      'nirK  –  copper nitrite reductase\n(NO₂ → NO)',           '#FFAB91'),
    ('norC/nosZ', 'norC / nosZ  –  NO & N₂O reductases\n(NO→N₂O→N₂)',       '#FFAB91'),
]
ROW_IDS = [r[0] for r in ROWS]

# ── Category colours ─────────────────────────────────────────────────────────
CAT_COLORS = {
    'GAO':                  '#E6811A',
    'PAO':                  '#2E86C1',
    'Denitrifier':          '#27AE60',
    'PAO/Denitrifier':      '#8E44AD',
    'GAO/Denitrifier':      '#D35400',
    'GAO/PAO':              '#F39C12',
    'GAO/PAO/Denitrifier':  '#7D3C98',
    'Other':                '#AAAAAA',
}
CAT_ORDER = ['GAO','PAO','Denitrifier','PAO/Denitrifier','GAO/Denitrifier',
             'GAO/PAO','GAO/PAO/Denitrifier','Other']

# ── 1. Load RPKM matrix ──────────────────────────────────────────────────────
rpkm = {}          # row_id -> {mag: value}
all_mags_ordered = []
with open(os.path.join(DATA, 'gene_rpkm.tsv')) as f:
    reader = csv.DictReader(f, delimiter='\t')
    all_mags_ordered = [c for c in reader.fieldnames if c != 'ko_group']
    for row in reader:
        rpkm[row['ko_group']] = {m: float(row[m]) for m in all_mags_ordered}

# ── 2. Load MAG metadata ─────────────────────────────────────────────────────
meta = {}   # mag -> {'classification','genus','species','iterid','abund_total'}
with open(os.path.join(DATA, 'selected_mags.tsv')) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mag = row['MAG']
        abund = sum(float(row[c]) for c in row if c.endswith('_relabund'))
        meta[mag] = {
            'classification': row['classification'],
            'abund': abund,
        }

# Taxonomy (loaded from taxonomy_labels.tsv if available)
tax_file = os.path.join(DATA, 'taxonomy_labels.tsv')
if os.path.exists(tax_file):
    with open(tax_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mag = row['MAG']
            if mag in meta:
                meta[mag]['genus'] = row.get('Genus', '')
                meta[mag]['species'] = row.get('Species', '')
                meta[mag]['iterid'] = row.get('mag_iterativeID', '')

# ── 3. Sort columns: by category then total abundance ───────────────────────
cat_priority = {c: i for i, c in enumerate(CAT_ORDER)}

def sort_key(mag):
    m = meta.get(mag, {})
    cat = m.get('classification', 'Other')
    p = cat_priority.get(cat, 99)
    return (p, -m.get('abund', 0))

sorted_mags = sorted(all_mags_ordered, key=sort_key)
n_cols = len(sorted_mags)
n_rows = len(ROWS)

# ── 4. Build matrix ──────────────────────────────────────────────────────────
matrix = np.zeros((n_rows, n_cols))
for i, (row_id, _, _) in enumerate(ROWS):
    for j, mag in enumerate(sorted_mags):
        matrix[i, j] = rpkm.get(row_id, {}).get(mag, 0)

log_matrix = np.log10(matrix + 1)
vmax = max(log_matrix.max(), 0.5)

# ── 5. Column labels ─────────────────────────────────────────────────────────
def fmt_label(mag):
    m = meta.get(mag, {})
    genus = m.get('genus', '')
    iterid = m.get('iterid', '')
    # Short MAG id: CAN_X_bin.N → CXbN, coasm_bin.N → cobN
    if mag.startswith('CAN_'):
        parts = mag.split('_')
        short = f"C{parts[1]}b{parts[-1].split('.')[-1]}"
    elif mag.startswith('coasm_bin.'):
        short = f"cob{mag.split('.')[-1]}"
    else:
        short = mag
    if genus:
        return f"{genus}\n({short})"
    elif iterid:
        base = iterid.split('.')[0].replace('Ca_','Ca.')
        return f"{base}\n({short})"
    return short

col_labels = [fmt_label(m) for m in sorted_mags]

# ── 6. Figure layout ─────────────────────────────────────────────────────────
col_w   = 0.5          # inches per column
row_h   = 0.65         # inches per row
lmargin = 4.0          # for row labels
rmargin = 1.0          # for colorbar
tmargin = 1.4          # for cat bar + legend
bmargin = 1.6          # for rotated col labels

fig_w = lmargin + n_cols * col_w + rmargin
fig_h = tmargin + n_rows * row_h + bmargin

fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

def frac(inches_from_left, width_inches, inches_from_bottom, height_inches):
    return [inches_from_left / fig_w, inches_from_bottom / fig_h,
            width_inches / fig_w, height_inches / fig_h]

heat_w  = n_cols * col_w
heat_h  = n_rows * row_h
heat_l  = lmargin
heat_b  = bmargin
cat_h   = 0.20
cat_b   = heat_b + heat_h + 0.04
strip_w = 0.12   # narrow colored strip between row labels and heatmap
strip_gap = 0.04

ax_heat  = fig.add_axes(frac(heat_l, heat_w, heat_b, heat_h))
ax_cat   = fig.add_axes(frac(heat_l, heat_w, cat_b,  cat_h))
ax_strip = fig.add_axes(frac(heat_l - strip_w - strip_gap, strip_w, heat_b, heat_h))

cmap = LinearSegmentedColormap.from_list(
    'rpkm', ['#FFFFFF', '#FEE8C8', '#FDBB84', '#E34A33', '#8B0000'], N=256)

# Heatmap
im = ax_heat.imshow(log_matrix, aspect='auto', cmap=cmap,
                    vmin=0, vmax=vmax, interpolation='nearest')

# Grid
for ax in [ax_heat]:
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which='minor', color='#E0E0E0', linewidth=0.5, zorder=2)
    ax.tick_params(which='minor', size=0)

# Row labels
ax_heat.set_yticks(range(n_rows))
ax_heat.set_yticklabels([r[1] for r in ROWS], fontsize=7, va='center')
ax_heat.tick_params(axis='y', length=0, pad=4)

# Row color strip (matches KO_groups.png color scheme)
for i, (_, _, color) in enumerate(ROWS):
    ax_strip.add_patch(mpatches.Rectangle((0, i - 0.5), 1, 1, color=color, linewidth=0))
ax_strip.set_xlim(0, 1)
ax_strip.set_ylim(-0.5, n_rows - 0.5)
ax_strip.set_xticks([])
ax_strip.set_yticks([])
for spine in ax_strip.spines.values():
    spine.set_visible(False)

# Column labels
ax_heat.set_xticks(range(n_cols))
ax_heat.set_xticklabels(col_labels, rotation=55, ha='right', fontsize=6.5, va='top')
ax_heat.tick_params(axis='x', length=0, pad=1)
ax_heat.set_xlabel('')

# Cell annotations (raw RPKM)
for i in range(n_rows):
    for j in range(n_cols):
        val = matrix[i, j]
        if val > 0:
            lv = log_matrix[i, j]
            tc = 'white' if lv > vmax * 0.65 else '#333333'
            txt = f'{val:.0f}' if val >= 10 else (f'{val:.1f}' if val >= 1 else f'{val:.2f}')
            ax_heat.text(j, i, txt, ha='center', va='center',
                         fontsize=4.5, color=tc, zorder=3)

# Category color bar
for j, mag in enumerate(sorted_mags):
    cat = meta.get(mag, {}).get('classification', 'Other')
    color = CAT_COLORS.get(cat, '#AAAAAA')
    ax_cat.add_patch(mpatches.FancyArrow(
        j + 0.5, 0.5, 0, 0, width=0.95, head_width=0, head_length=0, color=color))
    ax_cat.add_patch(mpatches.Rectangle((j + 0.025, 0.05), 0.95, 0.9,
                                         color=color, linewidth=0))
ax_cat.set_xlim(0, n_cols)
ax_cat.set_ylim(0, 1)
ax_cat.axis('off')

# Legend
seen_cats = [c for c in CAT_ORDER if any(
    meta.get(m, {}).get('classification', 'Other') == c for m in sorted_mags)]
patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in seen_cats]
ax_cat.legend(handles=patches, loc='lower center',
              bbox_to_anchor=(0.5, 1.05), ncol=len(patches),
              fontsize=7.5, frameon=False, handlelength=1.2)

# Relative abundance strip (thin bar below col labels — total sum across samples)
ax_ab = fig.add_axes(frac(heat_l, heat_w, bmargin - 0.28, 0.18))
abunds = np.array([meta.get(m, {}).get('abund', 0) for m in sorted_mags])
if abunds.max() > 0:
    abunds_norm = abunds / abunds.max()
    cmap_ab = plt.cm.YlOrBr
    for j, (a_norm, a_raw) in enumerate(zip(abunds_norm, abunds)):
        ax_ab.add_patch(mpatches.Rectangle(
            (j + 0.05, 0), 0.9, 1, color=cmap_ab(a_norm * 0.85 + 0.1)))
ax_ab.set_xlim(0, n_cols)
ax_ab.set_ylim(0, 1)
ax_ab.set_xticks([])
ax_ab.set_yticks([])
ax_ab.set_ylabel('Rel.abund.\n(sum)', fontsize=5, rotation=0, labelpad=28, va='center')
for spine in ax_ab.spines.values():
    spine.set_visible(False)

# RPKM colorbar
cb_ax = fig.add_axes(frac(heat_l + heat_w + 0.1, 0.15, heat_b + heat_h * 0.15,
                           heat_h * 0.70))
norm = Normalize(vmin=0, vmax=vmax)
cb = ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')
tick_log = np.linspace(0, vmax, 5)
cb.set_ticks(tick_log)
cb.set_ticklabels([f'{10**t - 1:.0f}' for t in tick_log])
cb_ax.set_ylabel('RPKM', fontsize=7)
cb_ax.tick_params(labelsize=6)

# Title
fig.text(0.5, 1.0 - 0.01,
         'Denitrification & DNRA Gene Abundance in Top High-Quality MAGs\n'
         '(long-read metagenome RPKM, all samples combined)',
         ha='center', va='top', fontsize=9, fontweight='bold',
         transform=fig.transFigure)

# Sample label
fig.text(heat_l / fig_w, (heat_b - 0.42) / fig_h,
         'CAN_1–5 samples (5 Nanopore metagenomes)',
         ha='left', va='bottom', fontsize=6.5, color='#666666',
         transform=fig.transFigure)

plt.savefig(os.path.join(WORK, 'gene_abundance_figure.png'),
            dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(WORK, 'gene_abundance_figure.pdf'),
            bbox_inches='tight', facecolor='white')
print(f"Saved: {WORK}/gene_abundance_figure.png")
print(f"Saved: {WORK}/gene_abundance_figure.pdf")
plt.close()
print("Done.")
