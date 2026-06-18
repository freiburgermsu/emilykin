#!/usr/bin/env python3
"""
Step 6: Plot gene copy number heatmap.

Mirrors 04_plot_heatmap.py layout:
  Rows  = KO functional groups
  Cols  = MAGs sorted by classification then total abundance
  Cells = gene copy number (linear scale; annotated with value)
  Top bar = classification colour
  Left strip = gene group colour

Copy number = (reads on gene / gene length) / (reads on genome / genome length).
A value of 1.0 means the gene is single-copy relative to genome coverage;
values > 1 indicate multiple copies or local enrichment.
"""

import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.colorbar import ColorbarBase
import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')

# ── Rows ─────────────────────────────────────────────────────────────────────
ROWS = [
    ('nirB/D',     'nirB/D  –  NADH nitrite reductase\n(DNRA: NO₂⁻ → NH₄⁺)',   '#A5D6A7'),
    ('nirS',       'nirS  –  cytochrome cd₁ nitrite reductase\n(denitr: NO₂⁻ → NO)',  '#A5D6A7'),
    ('nirK',       'nirK  –  copper nitrite reductase\n(denitr: NO₂⁻ → NO)',    '#FFAB91'),
    ('norBC/nosZ', 'norBC / nosZ  –  NO & N₂O reductases\n(denitr: NO→N₂O→N₂)', '#FFAB91'),
]
ROW_IDS = [r[0] for r in ROWS]

# ── Category colours ──────────────────────────────────────────────────────────
CAT_COLORS = {
    'GAO':                  '#E6811A',
    'PAO':                  '#2E86C1',
    'Denitrifier':          '#27AE60',
    'PAO/Denitrifier':      '#8E44AD',
    'Denitrifier/PAO':      '#8E44AD',
    'GAO/Denitrifier':      '#D35400',
    'GAO/PAO':              '#F39C12',
    'GAO/PAO/Denitrifier':  '#7D3C98',
    'Other':                '#AAAAAA',
}
CAT_ORDER = ['GAO', 'PAO', 'Denitrifier', 'PAO/Denitrifier', 'Denitrifier/PAO',
             'GAO/Denitrifier', 'GAO/PAO', 'GAO/PAO/Denitrifier', 'Other']

# ── 1. Load copy number matrix ────────────────────────────────────────────────
cn = {}               # row_id -> {mag: float}
all_mags_ordered = []
with open(os.path.join(DATA, 'gene_copy_number.tsv')) as f:
    reader = csv.DictReader(f, delimiter='\t')
    all_mags_ordered = [c for c in reader.fieldnames if c != 'ko_group']
    for row in reader:
        cn[row['ko_group']] = {m: float(row[m]) for m in all_mags_ordered}

# ── 2. Load MAG metadata ──────────────────────────────────────────────────────
meta = {}
with open(os.path.join(DATA, 'selected_mags.tsv')) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mag = row['MAG']
        abund = sum(float(row[c]) for c in row if c.endswith('_relabund'))
        meta[mag] = {'classification': row['classification'], 'abund': abund}

tax_file = os.path.join(DATA, 'taxonomy_labels.tsv')
if os.path.exists(tax_file):
    with open(tax_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mag = row['MAG']
            if mag in meta:
                meta[mag]['genus']   = row.get('Genus', '')
                meta[mag]['species'] = row.get('Species', '')
                meta[mag]['iterid']  = row.get('mag_iterativeID', '')

# ── 3. Sort columns: category priority, then total abundance ──────────────────
cat_priority = {c: i for i, c in enumerate(CAT_ORDER)}

def sort_key(mag):
    m   = meta.get(mag, {})
    cat = m.get('classification', 'Other')
    return (cat_priority.get(cat, 99), -m.get('abund', 0))

sorted_mags = sorted(all_mags_ordered, key=sort_key)
n_cols = len(sorted_mags)
n_rows = len(ROWS)

# ── 4. Build matrix ───────────────────────────────────────────────────────────
matrix = np.zeros((n_rows, n_cols))
for i, (row_id, _, _) in enumerate(ROWS):
    for j, mag in enumerate(sorted_mags):
        matrix[i, j] = cn.get(row_id, {}).get(mag, 0.0)

# Scale: log10(cn + 1) so that cn=0 → 0, cn=1 → 0.301, cn=10 → 1.041
log_matrix = np.log10(matrix + 1)
vmax = max(float(log_matrix.max()), 0.3)

# ── 5. Column labels ──────────────────────────────────────────────────────────
def fmt_label(mag):
    m      = meta.get(mag, {})
    genus  = m.get('genus', '')
    iterid = m.get('iterid', '')
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
        base = iterid.split('.')[0].replace('Ca_', 'Ca.')
        return f"{base}\n({short})"
    return short

col_labels = [fmt_label(m) for m in sorted_mags]

# ── 6. Figure layout (mirrors 04_plot_heatmap.py) ────────────────────────────
col_w   = 0.5
row_h   = 0.65
lmargin = 5.2
rmargin = 1.1
tmargin = 1.4
bmargin = 1.6

fig_w = lmargin + n_cols * col_w + rmargin
fig_h = tmargin + n_rows * row_h + bmargin

fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

def frac(l, w, b, h):
    return [l / fig_w, b / fig_h, w / fig_w, h / fig_h]

heat_w   = n_cols * col_w
heat_h   = n_rows * row_h
heat_l   = lmargin
heat_b   = bmargin
cat_h    = 0.20
cat_b    = heat_b + heat_h + 0.04
strip_w  = 0.12
strip_gap = 0.04

ax_heat  = fig.add_axes(frac(heat_l, heat_w, heat_b, heat_h))
ax_cat   = fig.add_axes(frac(heat_l, heat_w, cat_b, cat_h))
ax_strip = fig.add_axes(frac(heat_l - strip_w - strip_gap, strip_w, heat_b, heat_h))

# Blue-tinted colormap: white → steel-blue → navy
cmap = LinearSegmentedColormap.from_list(
    'cn', ['#FFFFFF', '#C6DBEF', '#6BAED6', '#2171B5', '#084594'], N=256)

im = ax_heat.imshow(log_matrix, aspect='auto', cmap=cmap,
                    vmin=0, vmax=vmax, interpolation='nearest')

# Grid
ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
ax_heat.grid(which='minor', color='#E0E0E0', linewidth=0.5, zorder=2)
ax_heat.tick_params(which='minor', size=0)

# Row labels — use figure-space Text so bbox_inches='tight' captures them correctly
ax_heat.set_yticks(range(n_rows))
ax_heat.set_yticklabels([], fontsize=7)
ax_heat.tick_params(axis='y', length=0, pad=0)
label_x_fig = (heat_l - strip_w - strip_gap - 0.08) / fig_w  # just left of colour strip
for i, (_, label, _) in enumerate(ROWS):
    # Convert data y-coord to figure fraction
    y_fig = (heat_b + (n_rows - 1 - i + 0.5) * row_h) / fig_h
    fig.text(label_x_fig, y_fig, label,
             ha='right', va='center', fontsize=7, transform=fig.transFigure)

# Row colour strip
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

# Cell annotations (raw copy number)
for i in range(n_rows):
    for j in range(n_cols):
        val = matrix[i, j]
        if val > 0.01:
            lv = log_matrix[i, j]
            tc = 'white' if lv > vmax * 0.60 else '#333333'
            txt = f'{val:.1f}'
            ax_heat.text(j, i, txt, ha='center', va='center',
                         fontsize=4.5, color=tc, zorder=3)

# Category colour bar
for j, mag in enumerate(sorted_mags):
    cat   = meta.get(mag, {}).get('classification', 'Other')
    color = CAT_COLORS.get(cat, '#AAAAAA')
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

# Relative abundance strip below column labels
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
ax_ab.set_ylabel('Rel. abund.\n(sum)', fontsize=5, rotation=0, labelpad=28, va='center')
for spine in ax_ab.spines.values():
    spine.set_visible(False)

# Colorbar — label in actual copy number units
cb_ax = fig.add_axes(frac(heat_l + heat_w + 0.12, 0.15,
                           heat_b + heat_h * 0.15, heat_h * 0.70))
norm = Normalize(vmin=0, vmax=vmax)
cb = ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')
# Tick positions in log space, labelled as real copy numbers
cn_ticks = [0, 0.5, 1, 2, 4, 8]
log_ticks = [np.log10(v + 1) for v in cn_ticks]
log_ticks_in_range = [(lt, v) for lt, v in zip(log_ticks, cn_ticks) if lt <= vmax * 1.01]
cb.set_ticks([lt for lt, _ in log_ticks_in_range])
cb.set_ticklabels([str(v) for _, v in log_ticks_in_range])
cb_ax.set_ylabel('Copy number\n(per genome)', fontsize=6.5, labelpad=2)
cb_ax.tick_params(labelsize=6)

# Reference line at cn=1 on colorbar
log_1 = np.log10(2)   # log10(1 + 1)
if log_1 <= vmax:
    cb_ax.axhline(log_1, color='#FF4444', linewidth=1.0, linestyle='--')
    cb_ax.text(1.08, log_1 / vmax, '1×', transform=cb_ax.transAxes,
               fontsize=5.5, color='#FF4444', va='center')

# Title
fig.text(0.5, 1.0 - 0.01,
         'Denitrification & DNRA Gene Copy Number in Top High-Quality MAGs\n'
         '(gene depth ÷ genome depth, all samples combined)',
         ha='center', va='top', fontsize=9, fontweight='bold',
         transform=fig.transFigure)

# Sample note
fig.text(heat_l / fig_w, (heat_b - 0.42) / fig_h,
         'CAN_1–5 samples (5 Nanopore metagenomes)',
         ha='left', va='bottom', fontsize=6.5, color='#666666',
         transform=fig.transFigure)

out_base = os.path.join(WORK, 'gene_copy_number_figure')
plt.savefig(out_base + '.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(out_base + '.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out_base}.png")
print(f"Saved: {out_base}.pdf")
print("Done.")
