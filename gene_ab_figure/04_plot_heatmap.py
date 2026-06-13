#!/usr/bin/env python3
"""
Gene abundance heatmap — denitrification / DNRA KOs across top HQ MAGs,
styled after example_gene_ab*.png.

  Columns = MAGs (union of top-10 HQ per sample), each MAG subdivided into the
            five Nanopore samples CAN_1–CAN_5 (left→right); boxed organism
            groups with italic genus headers.
  Rows    = KO functional groups, grouped by N-transformation stage with
            right-side brackets.
  Cells   = per-sample RPKM (log10 for colour); grey = no detected gene.
            KO-group cell per MAG/sample = max RPKM over the group's genes.
  Top bar = GAO / PAO / Denitrifier classification colour (qualitative palette,
            kept distinct from the coolwarm_r heatmap).

Data: data/gene_rpkm_per_sample.tsv, selected_mags.tsv, taxonomy_labels.tsv
Requires: matplotlib, numpy (+ stdlib csv). No seaborn.
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase
import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

# ── Row groups (N-transformation stage), top → bottom ───────────────────────
# (group label, [(row_id, row_label), ...])
ROW_GROUPS = [
    ('DNRA\n(NO₂→NH₄)',
     [('nirB/D', 'nirB/D  –  NADH nitrite reductase')]),
    ('Denitrification\n(NO₂→NO)',
     [('nirS', 'nirS  –  cytochrome cd₁ nitrite reductase'),
      ('nirK', 'nirK  –  copper nitrite reductase')]),
    ('Denitrification\n(NO→N₂O)',
     [('norC', 'norC  –  nitric oxide reductase')]),
    ('Denitrification\n(N₂O→N₂)',
     [('nosZ_cladeI',  'nosZ Clade I  –  N₂O reductase (typical)'),
      ('nosZ_cladeII', 'nosZ Clade II  –  N₂O reductase (atypical)')]),
]

# ── Column category colours (qualitative; distinct from coolwarm_r) ─────────
CAT_COLORS = {
    'PAO':                 '#7B3294',   # purple
    'Denitrifier/PAO':     '#1B9E77',   # teal-green
    'Denitrifier':         '#E6AB02',   # gold
    'GAO':                 '#D95F02',   # orange
    'GAO/PAO':             '#E7298A',   # magenta
    'Denitrifier/GAO':     '#A6761D',   # brown
    'Denitrifier/GAO/PAO': '#66A61E',   # olive
    'Other':               '#999999',   # grey
}
_CANON = ['Denitrifier', 'GAO', 'PAO']
def canon_cat(cat):
    cat = (cat or 'Other').strip()
    if cat in CAT_COLORS:
        return cat
    parts = [p.strip() for p in cat.replace('+', '/').split('/') if p.strip()]
    if not parts:
        return 'Other'
    parts = sorted(set(parts), key=lambda p: _CANON.index(p) if p in _CANON else 9)
    key = '/'.join(parts)
    return key if key in CAT_COLORS else 'Other'
CAT_ORDER = ['PAO', 'Denitrifier/PAO', 'Denitrifier', 'GAO/PAO',
             'Denitrifier/GAO', 'Denitrifier/GAO/PAO', 'GAO', 'Other']

# ── Load per-sample RPKM ────────────────────────────────────────────────────
psm, mags = {}, []
with open(os.path.join(DATA, 'gene_rpkm_per_sample_claded.tsv')) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g, m = row['ko_group'], row['mag']
        psm.setdefault(g, {}).setdefault(m, {})
        for s in SAMPLES:
            psm[g][m][s] = float(row[s])
        if m not in mags:
            mags.append(m)

# ── MAG metadata ────────────────────────────────────────────────────────────
meta = {}
with open(os.path.join(DATA, 'selected_mags.tsv')) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        m = row['MAG']
        meta[m] = {'classification': canon_cat(row['classification']),
                   'abund': sum(float(row[c]) for c in row if c.endswith('_relabund'))}
taxf = os.path.join(DATA, 'taxonomy_labels.tsv')
if os.path.exists(taxf):
    with open(taxf) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['MAG'] in meta:
                meta[row['MAG']]['genus'] = row.get('Genus', '')
                meta[row['MAG']]['iterid'] = row.get('mag_iterativeID', '')

# nosZ clade calls (for low-confidence markers)
nosz_clade, nosz_lowconf = {}, set()
czp = os.path.join(DATA, 'nosz_clades.tsv')
if os.path.exists(czp):
    with open(czp) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            nosz_clade[r['mag']] = r['clade']
            if r.get('confidence') == 'low':
                nosz_lowconf.add(r['mag'])

cat_priority = {c: i for i, c in enumerate(CAT_ORDER)}
def sort_key(m):
    mm = meta.get(m, {})
    return (cat_priority.get(mm.get('classification', 'Other'), 99), -mm.get('abund', 0))
sorted_mags = sorted(mags, key=sort_key)
n_mag, n_s = len(sorted_mags), len(SAMPLES)

def short_id(m):
    if m.startswith('CAN_'):
        p = m.split('_'); return f"C{p[1]}b{p[-1].split('.')[-1]}"
    if m.startswith('coasm_bin.'):
        return f"cob{m.split('.')[-1]}"
    return m
def genus_of(m):
    g = meta.get(m, {}).get('genus', '')
    if g:
        return g
    it = meta.get(m, {}).get('iterid', '')
    return it.split('.')[0].replace('Ca_', 'Ca.') if it else short_id(m)

# ── Geometry (data units) ───────────────────────────────────────────────────
SUB_W, GAP_X = 1.0, 0.6
ROW_H, GAP_Y = 1.0, 0.55
mag_x0, x = {}, 0.0
for m in sorted_mags:
    mag_x0[m] = x; x += n_s * SUB_W + GAP_X
total_x = x - GAP_X

row_y, row_labels, group_spans, y = {}, {}, [], 0.0
for glabel, rows in ROW_GROUPS:
    ytop = y
    for rid, rlab in rows:
        row_y[rid] = y; row_labels[rid] = rlab; y += ROW_H
    group_spans.append((glabel, ytop, y))
    y += GAP_Y
total_y = y - GAP_Y

# ── Colour scale: coolwarm_r over present values; empty = light grey ────────
allv = [psm[g][m][s] for g in psm for m in psm[g] for s in SAMPLES]
vmax = np.log10(max(allv) + 1)
norm = Normalize(vmin=0, vmax=vmax)
# light → dark blue sequential map, matching this project's relative-abundance
# heatmap (data_processing.py / render_abundance_heatmaps.py)
cmap = LinearSegmentedColormap.from_list(
    'NewMap', [(0.0, 'aliceblue'), (0.25, 'lightblue'), (1.0, 'navy')])
EMPTY = '#E8E8E8'

# ── Figure ───────────────────────────────────────────────────────────────────
xs, ys = 0.145, 0.62                      # inches per data unit
lmargin, rmargin = 4.3, 2.9
tmargin, bmargin = 1.3, 2.7
plot_w, plot_h = total_x * xs, total_y * ys
fig_w = lmargin + plot_w + rmargin
fig_h = tmargin + plot_h + bmargin
fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

def frac(l, w, b, h):
    return [l / fig_w, b / fig_h, w / fig_w, h / fig_h]

ax = fig.add_axes(frac(lmargin, plot_w, bmargin, plot_h))
ax.set_xlim(0, total_x); ax.set_ylim(total_y, 0)      # invert y → row 0 at top
ax.axis('off')

# cells
for rid, y0 in row_y.items():
    for m in sorted_mags:
        for si, s in enumerate(SAMPLES):
            v = psm[rid][m][s]
            c = cmap(norm(np.log10(v + 1))) if v > 0 else EMPTY
            ax.add_patch(Rectangle((mag_x0[m] + si * SUB_W, y0), SUB_W, ROW_H,
                                   facecolor=c, edgecolor='white', linewidth=0.4))

# boxed organism groups (full height)
for m in sorted_mags:
    ax.add_patch(Rectangle((mag_x0[m], 0), n_s * SUB_W, total_y, fill=False,
                           edgecolor='#555555', linewidth=0.8, zorder=4))

# low-confidence nosZ clade-call markers (hollow circle on the clade row)
for m in nosz_lowconf:
    rid = 'nosZ_cladeI' if nosz_clade.get(m) == 'I' else 'nosZ_cladeII'
    if rid in row_y and m in mag_x0:
        ax.plot(mag_x0[m] + n_s * SUB_W / 2, row_y[rid] + ROW_H / 2, marker='o',
                mfc='none', mec='#111111', mew=1.0, ms=6, clip_on=False, zorder=6)

# row labels (left)
for rid, y0 in row_y.items():
    ax.text(-0.4, y0 + ROW_H / 2, row_labels[rid], ha='right', va='center',
            fontsize=8, clip_on=False)

# right-side N-stage brackets
xb = total_x + 0.4
for glabel, ytop, ybot in group_spans:
    ax.plot([xb, xb + 0.45, xb + 0.45, xb],
            [ytop + 0.06, ytop + 0.06, ybot - 0.06, ybot - 0.06],
            color='#444444', lw=1.1, clip_on=False)
    ax.text(xb + 0.75, (ytop + ybot) / 2, glabel, ha='left', va='center',
            rotation=90, fontsize=8, clip_on=False)

# bottom: per-sub-column sample index, per-block italic genus + short id
for m in sorted_mags:
    bx = mag_x0[m]
    for si in range(n_s):
        ax.text(bx + si * SUB_W + 0.5, total_y + 0.18, str(si + 1),
                ha='center', va='top', fontsize=4.6, color='#555555', clip_on=False)
    ax.text(bx + n_s * SUB_W / 2, total_y + 0.65,
            f"{genus_of(m)}\n({short_id(m)})", ha='right', va='top',
            rotation=45, rotation_mode='anchor', fontsize=6.6,
            fontstyle='italic', clip_on=False)

# ── Category colour bar + legend (top) ──────────────────────────────────────
catb = bmargin + plot_h + 0.05
axc = fig.add_axes(frac(lmargin, plot_w, catb, 0.17))
axc.set_xlim(0, total_x); axc.set_ylim(0, 1); axc.axis('off')
for m in sorted_mags:
    col = CAT_COLORS.get(meta.get(m, {}).get('classification', 'Other'), '#999999')
    axc.add_patch(Rectangle((mag_x0[m], 0.1), n_s * SUB_W, 0.8,
                            facecolor=col, edgecolor='none'))
seen = [c for c in CAT_ORDER
        if any(meta.get(m, {}).get('classification', 'Other') == c for m in sorted_mags)]
handles = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in seen]
handles.append(mpatches.Patch(color=EMPTY, label='no gene detected'))
axc.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.7),
           ncol=len(handles), fontsize=8, frameon=False, handlelength=1.3,
           columnspacing=1.4)

# ── RPKM colorbar (far right) ───────────────────────────────────────────────
cb_ax = fig.add_axes(frac(lmargin + plot_w + rmargin - 0.55, 0.16,
                          bmargin + plot_h * 0.2, plot_h * 0.6))
cb = ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')
ticks = np.linspace(0, vmax, 5)
cb.set_ticks(ticks)
cb.set_ticklabels([f'{10**t - 1:.0f}' for t in ticks])
cb_ax.set_ylabel('RPKM', fontsize=8)
cb_ax.tick_params(labelsize=6.5)

# ── Title + serif caption ───────────────────────────────────────────────────
fig.text(lmargin / fig_w, 1 - 0.015,
         'Denitrification & DNRA gene abundance across top high-quality MAGs '
         'of the CAN system',
         ha='left', va='top', fontsize=11, fontweight='bold')

caption = (
    'Per-sample gene abundances in reads per kilobase per million mapped reads '
    '(RPKM) from long-read (Nanopore) metagenomes of the CAN system. Columns are '
    'grouped by MAG — the union of the ten most abundant high-quality MAGs per '
    'sample — each subdivided into the five samples CAN_1–CAN_5 (sub-columns 1–5, '
    'left→right); the coloured bar marks GAO / PAO / denitrifier classification '
    '(mag_abundance_summary.tsv logic). Rows are denitrification/DNRA KO functions '
    'grouped by nitrogen-transformation stage (right brackets); for nirB/D '
    '(K00362/K00363) the cell is the maximum RPKM over the two genes. nosZ (K00376) '
    'is split into Clade I (typical denitrifier) and Clade II (atypical / '
    'non-denitrifier) by per-gene nhmmer classification against the NosZREF Clade I '
    'and Clade II (sub-clades 1577 A–H) reference alignments; ○ marks a '
    'low-confidence clade call. Grey cells indicate no detected gene.')
fig.text(lmargin / fig_w, (bmargin - 1.05) / fig_h, caption,
         ha='left', va='top', fontsize=7.6, family='serif', wrap=True)
# constrain caption wrap width
fig.texts[-1]._get_wrap_line_width = lambda: (plot_w + rmargin) * fig.dpi

plt.savefig(os.path.join(WORK, 'gene_abundance_figure.png'),
            dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(WORK, 'gene_abundance_figure.pdf'),
            bbox_inches='tight', facecolor='white')
print(f"Saved gene_abundance_figure.png / .pdf  "
      f"({n_mag} MAGs × {n_s} samples = {n_mag*n_s} sub-columns, "
      f"vmax={10**vmax-1:.0f} RPKM)")
plt.close()
print("Done.")
