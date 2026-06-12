#!/usr/bin/env python3
"""Quick test: generate figure with proxy RPKM from coverm short-read depths.
Useful to verify layout before long-read alignment completes.
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.colorbar import ColorbarBase
import numpy as np

WORK = '/scratch1/afreiburger/emilykin/gene_ab_figure'
DATA = f'{WORK}/data'
PROC = '/scratch1/afreiburger/emilykin/processed/mag'
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

KO_GROUPS = [
    ('nirB/D',    ['K00362', 'K00363']),
    ('nirS',      ['K15864']),
    ('nirK',      ['K00368']),
    ('norC/nosZ', ['K02305', 'K00376']),
]

# ── Load MAG metadata ────────────────────────────────────────────────────────
selected_mags = [l.strip() for l in open(f'{DATA}/selected_mag_list.txt') if l.strip()]

meta = {}
with open(f'{DATA}/selected_mags.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        meta[row['MAG']] = {'classification': row['classification'],
                            'abund': sum(float(row[c]) for c in row if c.endswith('_relabund'))}

# Load taxonomy
with open(f'{DATA}/taxonomy_labels.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['MAG'] in meta:
            meta[row['MAG']].update({'genus': row['Genus'],
                                     'species': row['Species'],
                                     'iterid': row['mag_iterativeID']})

# ── Load coverm trimmed mean depths (proxy for gene presence + level) ────────
coverm_depth = {}  # mag -> sample -> trimmed_mean
with open(f'{PROC}/abundance/coverm_genome.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mag = row['Genome']
        if mag == 'unmapped': continue
        coverm_depth[mag] = {s: float(row.get(f'{s} Trimmed Mean', 0) or 0)
                             for s in SAMPLES}

# ── Load target genes (which KO belongs to which MAG) ───────────────────────
mag_kos = {m: set() for m in selected_mags}
with open(f'{DATA}/target_genes.bed') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mag_kos[row['mag']].add(row['ko'])

# ── Build proxy RPKM matrix ──────────────────────────────────────────────────
# Proxy: if gene is present in MAG, use sum of trimmed mean depths across samples
# If gene not detected, 0.
ROWS = [
    ('nirB/D',    'nirB/D  –  NADH nitrite reductase\n(DNRA: NO₂ → NH₄)'),
    ('nirS',      'nirS  –  cytochrome cd₁ nitrite reductase\n(NO₂ → NO)'),
    ('nirK',      'nirK  –  copper nitrite reductase\n(NO₂ → NO)'),
    ('norC/nosZ', 'norC / nosZ  –  NO & N₂O reductases\n(NO→N₂O→N₂)'),
]
ROW_IDS = [r[0] for r in ROWS]

KO_TO_GROUP = {}
for gname, kos in KO_GROUPS:
    for ko in kos:
        KO_TO_GROUP[ko] = gname

matrix = {}  # group -> mag -> value
for gname, kos in KO_GROUPS:
    matrix[gname] = {}
    for mag in selected_mags:
        mag_kos_set = mag_kos.get(mag, set())
        if any(ko in mag_kos_set for ko in kos):
            # Gene present: use sum of trimmed mean across samples as proxy
            depth_sum = sum(coverm_depth.get(mag, {}).get(s, 0) for s in SAMPLES)
            matrix[gname][mag] = depth_sum
        else:
            matrix[gname][mag] = 0.0

# ── Sort columns ─────────────────────────────────────────────────────────────
CAT_COLORS = {
    'GAO': '#E6811A', 'PAO': '#2E86C1', 'Denitrifier': '#27AE60',
    'PAO/Denitrifier': '#8E44AD', 'GAO/Denitrifier': '#D35400',
    'GAO/PAO': '#F39C12', 'GAO/PAO/Denitrifier': '#7D3C98', 'Other': '#AAAAAA',
}
CAT_ORDER = ['GAO', 'PAO', 'Denitrifier', 'PAO/Denitrifier', 'GAO/Denitrifier',
             'GAO/PAO', 'GAO/PAO/Denitrifier', 'Other']
cat_pri = {c: i for i, c in enumerate(CAT_ORDER)}

def sort_key(mag):
    m = meta.get(mag, {})
    return (cat_pri.get(m.get('classification', 'Other'), 99), -m.get('abund', 0))

sorted_mags = sorted(selected_mags, key=sort_key)
n_cols = len(sorted_mags)
n_rows = len(ROWS)

mat = np.zeros((n_rows, n_cols))
for i, (rid, _) in enumerate(ROWS):
    for j, mag in enumerate(sorted_mags):
        mat[i, j] = matrix[rid].get(mag, 0)

log_mat = np.log10(mat + 1)
vmax = max(log_mat.max(), 0.5)

# Column labels
def fmt_label(mag):
    m = meta.get(mag, {})
    genus = m.get('genus', '')
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

# ── Plot ─────────────────────────────────────────────────────────────────────
col_w, row_h = 0.5, 0.65
lmargin, rmargin, tmargin, bmargin = 4.0, 1.0, 1.4, 1.8

fig_w = lmargin + n_cols * col_w + rmargin
fig_h = tmargin + n_rows * row_h + bmargin
fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

def frac(l, w, b, h):
    return [l/fig_w, b/fig_h, w/fig_w, h/fig_h]

heat_w, heat_h = n_cols * col_w, n_rows * row_h
cat_h = 0.22

ax_heat = fig.add_axes(frac(lmargin, heat_w, bmargin, heat_h))
ax_cat  = fig.add_axes(frac(lmargin, heat_w, bmargin + heat_h + 0.04, cat_h))

cmap = LinearSegmentedColormap.from_list(
    'rpkm', ['#FFFFFF', '#FEE8C8', '#FDBB84', '#E34A33', '#8B0000'], N=256)

im = ax_heat.imshow(log_mat, aspect='auto', cmap=cmap, vmin=0, vmax=vmax,
                    interpolation='nearest')

# Grid
ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
ax_heat.grid(which='minor', color='#E0E0E0', linewidth=0.5, zorder=2)
ax_heat.tick_params(which='minor', size=0)

ax_heat.set_yticks(range(n_rows))
ax_heat.set_yticklabels([r[1] for r in ROWS], fontsize=7.5, va='center')
ax_heat.tick_params(axis='y', length=0, pad=4)

ax_heat.set_xticks(range(n_cols))
ax_heat.set_xticklabels(col_labels, rotation=55, ha='right', fontsize=7, va='top')
ax_heat.tick_params(axis='x', length=0, pad=1)

# Cell values
for i in range(n_rows):
    for j in range(n_cols):
        val = mat[i, j]
        if val > 0:
            lv = log_mat[i, j]
            tc = 'white' if lv > vmax * 0.65 else '#333333'
            txt = f'{val:.0f}'
            ax_heat.text(j, i, txt, ha='center', va='center',
                         fontsize=4.0, color=tc, zorder=3)

# Category bar
for j, mag in enumerate(sorted_mags):
    cat = meta.get(mag, {}).get('classification', 'Other')
    color = CAT_COLORS.get(cat, '#AAAAAA')
    ax_cat.add_patch(mpatches.Rectangle((j + 0.02, 0.05), 0.96, 0.9, color=color))
ax_cat.set_xlim(0, n_cols)
ax_cat.set_ylim(0, 1)
ax_cat.axis('off')

seen_cats = [c for c in CAT_ORDER if any(
    meta.get(m, {}).get('classification', 'Other') == c for m in sorted_mags)]
patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in seen_cats]
ax_cat.legend(handles=patches, loc='lower center',
              bbox_to_anchor=(0.5, 1.05), ncol=len(patches),
              fontsize=7.5, frameon=False, handlelength=1.2)

# Relative abundance strip
ax_ab = fig.add_axes(frac(lmargin, heat_w, bmargin - 0.3, 0.2))
abunds = np.array([meta.get(m, {}).get('abund', 0) for m in sorted_mags])
if abunds.max() > 0:
    abunds_n = abunds / abunds.max()
    cmap_ab = plt.cm.YlOrBr
    for j, a in enumerate(abunds_n):
        ax_ab.add_patch(mpatches.Rectangle((j + 0.05, 0), 0.9, 1,
                                            color=cmap_ab(a * 0.85 + 0.1)))
ax_ab.set_xlim(0, n_cols)
ax_ab.set_ylim(0, 1)
ax_ab.axis('off')
ax_ab.set_ylabel('Rel.abund', fontsize=5.5, rotation=0, labelpad=30, va='center')

# Colorbar
cb_ax = fig.add_axes(frac(lmargin + heat_w + 0.12, 0.15,
                           bmargin + heat_h * 0.15, heat_h * 0.7))
norm = Normalize(vmin=0, vmax=vmax)
cb = ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')
ticks = np.linspace(0, vmax, 5)
cb.set_ticks(ticks)
cb.set_ticklabels([f'{10**t - 1:.0f}' for t in ticks])
cb_ax.set_ylabel('Proxy value\n(short-read\ntrimmed mean\ndepth, sum)', fontsize=5.5)
cb_ax.tick_params(labelsize=5.5)

fig.text(0.5, 1.0 - 0.005,
         'Denitrification & DNRA Gene Abundance in Top HQ MAGs\n'
         '(PROXY: short-read trimmed mean depth — waiting for long-read alignment)',
         ha='center', va='top', fontsize=8.5, fontweight='bold',
         transform=fig.transFigure, color='#444444')

out = f'{WORK}/gene_abundance_figure_PROXY.png'
plt.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}")
plt.close()
