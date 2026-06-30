#!/usr/bin/env python3
"""
Denitrification gene inventory of the NOVEL bacteria — Figure-S3-style bubble plot.

Reference: Figure S3 of the Accumulibacter metatranscriptomics paper (image.png)
plots RNA-RPKM gene *expression* across redox conditions.  This project has no
metatranscriptome; its denitrification data are *metagenomic*.  So this is the
faithful metagenomic analog that DEFINES each novel MAG's denitrification
potential:

  rows    = denitrification genes, grouped into the same numbered N-reduction
            steps as Figure S3 (1 nitrate, 2 nitrite, 3 nitric-oxide, 4 nitrous-
            oxide; plus the DNRA branch).  Step "0 - nitrate/nitrite transporter"
            is omitted: no narK/nrt transporter gene was among the target loci.
  columns = the species-novel focal MAGs (no GTDB species assignment / newly
            named Candidatus), labelled by iterativeID, grouped/headed by genus,
            with a phylum sub-label and an N-cycle classification colour strip.
  bubble SIZE  = mean DNA-RPKM of the gene group across the five metagenomes
                 CAN_1-CAN_5 (community gene-copy abundance; the analog of the
                 reference's transcript abundance, averaged across cycles).
  bubble COLOUR= gene copies per genome (the genetic "definition" of the
                 organism's denitrification potential).
  faint dot    = gene absent (copy number 0).

Data (all in ./data, produced by the gene_ab_figure pipeline):
  gene_rpkm_per_sample_claded.tsv   per-group RPKM x MAG x {CAN_1..5}
  gene_copy_per_mag.tsv             per-gene copy number + iterativeID + genus
  selected_mags.tsv                 N-cycle classification (Denitrifier/PAO/...)
  taxonomy_labels.tsv               GTDB Genus/Species (novelty = blank Species)
  ../gtdbtk.bac120.summary.tsv      GTDB phylum/class for the genus sub-label
Requires: matplotlib, numpy (+ stdlib csv).  No seaborn.
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import BoundaryNorm
from matplotlib.lines import Line2D
import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')
ROOT = os.path.dirname(WORK)
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

# ── Figure rows: (group_id_in_RPKM, gene_label, copy-number member columns) ──
# copy number for a group = max over its subunit columns (operon copy count).
ROWS = [
    ('nar',          'nar  (narGHI)',  ['narG', 'narH', 'narI']),
    ('nap',          'nap  (napAB)',   ['napA', 'napB']),
    ('nirS',         'nirS',           ['nirS']),
    ('nirK',         'nirK',           ['nirK']),
    ('nirB/D',       'nirB/D',         ['nirB', 'nirD']),
    ('norC',         'nor  (norBC)',   ['norB', 'norC']),
    ('nosZ_cladeI',  'nosZ  clade I',  ['nosZ_cladeI']),
    ('nosZ_cladeII', 'nosZ  clade II', ['nosZ_cladeII']),
]
ROW_IDS = [r[0] for r in ROWS]

# Right-side step brackets (Figure-S3 numbering): (label, first_row, last_row)
STEPS = [
    ('1 · Nitrate reduction\n(NO$_3^-\\!\\to$NO$_2^-$)',      0, 1),
    ('2 · Nitrite reduction\n(NO$_2^-\\!\\to$NO)',           2, 3),
    ('DNRA\n(NO$_2^-\\!\\to$NH$_4^+$)',                       4, 4),
    ('3 · Nitric-oxide reduction\n(NO$\\to$N$_2$O)',          5, 5),
    ('4 · Nitrous-oxide reduction\n(N$_2$O$\\to$N$_2$)',      6, 7),
]

# ── N-cycle classification colours (match 04_plot_heatmap.py palette) ───────
CAT_COLORS = {
    'PAO': '#7B3294', 'Denitrifier/PAO': '#1B9E77', 'Denitrifier': '#E6AB02',
    'GAO': '#D95F02', 'GAO/PAO': '#E7298A', 'Denitrifier/GAO': '#A6761D',
    'Denitrifier/GAO/PAO': '#66A61E', 'Other': '#999999',
}

def load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

# ── 1. RPKM: (group, mag) -> mean over the 5 metagenomes ────────────────────
rpkm_mean = {}
for r in load_tsv(os.path.join(DATA, 'gene_rpkm_per_sample_claded.tsv')):
    vals = [float(r[s]) for s in SAMPLES]
    rpkm_mean[(r['ko_group'], r['mag'])] = float(np.mean(vals))

# ── 2. Copy number + iterativeID + genus ────────────────────────────────────
copy, iterid, genus = {}, {}, {}
for r in load_tsv(os.path.join(DATA, 'gene_copy_per_mag.tsv')):
    m = r['MAG']
    iterid[m] = r['iterativeID']
    genus[m] = r['genus']
    copy[m] = {k: int(float(r[k])) for k in r
               if k not in ('MAG', 'iterativeID', 'genus', 'total',
                            'nosZ_II_subclade', 'ppk1_type', 'ppk1_species')}

# ── 3. classification + GTDB species (novelty) ──────────────────────────────
classification, species = {}, {}
for r in load_tsv(os.path.join(DATA, 'selected_mags.tsv')):
    classification[r['MAG']] = r['classification']
for r in load_tsv(os.path.join(DATA, 'taxonomy_labels.tsv')):
    species[r['MAG']] = r['Species'].strip()

# ── 4. GTDB phylum/class for the genus sub-label ────────────────────────────
phylum = {}
for r in load_tsv(os.path.join(ROOT, 'gtdbtk.bac120.summary.tsv')):
    d = dict(x.split('__', 1) for x in r['classification'].split(';'))
    p = d.get('p', '').strip()
    p = {'Pseudomonadota': 'Proteobacteria'}.get(p, p)  # GTDB->common
    phylum[r['user_genome']] = p

# ── 5. Define the NOVEL set: focal MAGs with no GTDB species assignment ──────
# (`classification` only holds the 21 real focal MAGs, so this drops the TOTAL row)
focal = [m for m in copy if m in classification]
novel = [m for m in focal if not species.get(m, '')]

def group_copy(m, members):
    return max((copy[m].get(c, 0) for c in members), default=0)

def total_denit(m):
    return sum(group_copy(m, mem) for _, _, mem in ROWS)

# columns: most-complete denitrifier (left) -> N2O-consumer only (right)
novel.sort(key=lambda m: (-total_denit(m), iterid[m]))

print(f'Focal MAGs: {len(focal)};  species-novel (blank GTDB species): {len(novel)}')
for m in novel:
    print(f'  {iterid[m]:26s} {genus[m]:18s} {phylum.get(m,"?"):16s} '
          f'tot_denit_copies={total_denit(m)}  class={classification.get(m,"?")}')

# ── Build the matrix ────────────────────────────────────────────────────────
ncol, nrow = len(novel), len(ROWS)
size = np.full((nrow, ncol), np.nan)   # mean RPKM where gene present
cval = np.full((nrow, ncol), np.nan)   # copy number where gene present
for j, m in enumerate(novel):
    for i, (gid, _, members) in enumerate(ROWS):
        cn = group_copy(m, members)
        if cn > 0:
            cval[i, j] = cn
            size[i, j] = rpkm_mean.get((gid, m), 0.0)

# ── Plot ────────────────────────────────────────────────────────────────────
plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(13.5, 8.4))

# faint background dot at every cell (present or not), like Figure S3
for i in range(nrow):
    for j in range(ncol):
        ax.scatter(j, i, s=7, c='#d9d9d9', zorder=1, edgecolors='none')

# bubble size scaling: area proportional to mean RPKM
SCALE = 13.0           # pts^2 per RPKM unit
MINSZ = 22.0           # floor so a present-but-near-zero gene stays visible
def bubble_s(v):
    return max(v * SCALE, MINSZ)

# copy-number colour: discrete Blues, 1/2/3+
cmap = plt.get_cmap('Blues')
maxcopy = int(np.nanmax(cval)) if np.isfinite(np.nanmax(cval)) else 1
bounds = [0.5, 1.5, 2.5, max(3.5, maxcopy + 0.5)]
norm = BoundaryNorm(bounds, cmap.N)
disc = {1: cmap(0.40), 2: cmap(0.70), 3: cmap(0.95)}

for i in range(nrow):
    for j in range(ncol):
        if np.isfinite(cval[i, j]):
            cn = int(cval[i, j])
            ax.scatter(j, i, s=bubble_s(size[i, j]),
                       color=disc.get(min(cn, 3)),
                       edgecolors='#2c3e50', linewidths=0.7, zorder=3)

# axes framing -------------------------------------------------------------
ax.set_xlim(-0.6, ncol - 0.4)
ax.set_ylim(-0.6, nrow - 0.5)
ax.invert_yaxis()                      # row 0 (nar) at top
ax.set_xticks(range(ncol))
ax.set_xticklabels([iterid[m] for m in novel], rotation=40, ha='right', fontsize=9)
ax.set_yticks(range(nrow))
ax.set_yticklabels([r[1] for r in ROWS], fontsize=9.5, style='italic')
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
# light cell gridlines
for x in np.arange(-0.5, ncol, 1):
    ax.axvline(x, color='#eeeeee', lw=0.8, zorder=0)
for y in np.arange(-0.5, nrow, 1):
    ax.axhline(y, color='#eeeeee', lw=0.8, zorder=0)

# right-side step brackets --------------------------------------------------
xb = ncol - 0.30
for label, r0, r1 in STEPS:
    ax.plot([xb, xb], [r0 - 0.32, r1 + 0.32], color='#444444', lw=1.3,
            clip_on=False, zorder=5)
    for yy in (r0 - 0.32, r1 + 0.32):
        ax.plot([xb - 0.06, xb], [yy, yy], color='#444444', lw=1.3,
                clip_on=False, zorder=5)
    ax.text(xb + 0.14, (r0 + r1) / 2, label, va='center', ha='left',
            fontsize=8.3, clip_on=False)

# top: phylum headers + N-cycle classification colour strip -----------------
# (genus is already carried by the iterativeID on the x-axis, so the top only
#  adds taxonomic context: the novel N-cyclers span four phyla.)
strip_y = -0.5
for j, m in enumerate(novel):
    ax.add_patch(Rectangle((j - 0.42, strip_y - 0.28), 0.84, 0.24,
                           color=CAT_COLORS.get(classification.get(m, 'Other'),
                                                '#999999'),
                           clip_on=False, zorder=4))
    ax.text(j, strip_y - 0.34, phylum.get(m, ''), rotation=45, ha='left',
            va='bottom', fontsize=8.0, color='#666666', clip_on=False)

# ── Legends (size + copy-number colour + classification) ────────────────────
# bubble-size legend
size_vals = [10, 25, 50]
size_handles = [Line2D([], [], marker='o', linestyle='', markersize=np.sqrt(bubble_s(v)),
                       markerfacecolor='#bdbdbd', markeredgecolor='#2c3e50',
                       markeredgewidth=0.7) for v in size_vals]
leg1 = ax.legend(size_handles, [f'{v}' for v in size_vals],
                 title='mean DNA-RPKM\n(gene abundance)', loc='upper left',
                 bbox_to_anchor=(1.36, 1.00), labelspacing=1.6, frameon=False,
                 handletextpad=1.1, borderpad=0.8, fontsize=8.5, title_fontsize=8.8)
ax.add_artist(leg1)

# copy-number colour legend
copy_handles = [Line2D([], [], marker='o', linestyle='', markersize=9,
                       markerfacecolor=disc[c], markeredgecolor='#2c3e50',
                       markeredgewidth=0.7) for c in (1, 2, 3)]
leg2 = ax.legend(copy_handles, ['1', '2', '3+'],
                 title='gene copies\nper genome', loc='upper left',
                 bbox_to_anchor=(1.36, 0.52), labelspacing=0.7, frameon=False,
                 fontsize=8.5, title_fontsize=8.8)
ax.add_artist(leg2)

# classification legend (only categories present)
present_cats = []
for m in novel:
    c = classification.get(m, 'Other')
    if c not in present_cats:
        present_cats.append(c)
cat_handles = [Rectangle((0, 0), 1, 1, color=CAT_COLORS.get(c, '#999999'))
               for c in present_cats]
leg3 = ax.legend(cat_handles, present_cats, title='N-cycle role', loc='upper left',
                 bbox_to_anchor=(1.36, 0.22), frameon=False, fontsize=8.3,
                 title_fontsize=8.8, handlelength=1.1, handleheight=1.1)
ax.add_artist(leg3)

# title + caption -----------------------------------------------------------
fig.suptitle('Denitrification gene inventory of the novel bacteria',
             x=0.035, y=0.975, ha='left', fontsize=14, fontweight='bold')
fig.text(0.035, 0.935,
         'species-novel MAGs (no GTDB species assignment) — metagenomic analog of Figure S3',
         ha='left', fontsize=9.5, color='#555555')

cap = ('Bubble area = mean DNA-RPKM of the gene group across the five metagenomes '
       'CAN_1–CAN_5 (community gene-copy abundance; analog of the reference\'s '
       'transcript abundance, averaged across samples). Colour = gene copies per '
       'genome. Faint grey dots mark absent genes. Steps follow Figure S3 '
       '(1 nitrate, 2 nitrite, 3 nitric-oxide, 4 nitrous-oxide reduction; DNRA '
       'branch shown separately); the "0 transporter" step is omitted — no '
       'narK/nrt gene was among the target loci. This figure is metagenomic gene '
       'abundance, NOT transcriptomic expression.')
fig.text(0.035, 0.015, cap, ha='left', va='bottom', fontsize=7.2,
         color='#555555', wrap=True)

fig.subplots_adjust(left=0.095, right=0.58, top=0.80, bottom=0.235)
out_png = os.path.join(WORK, 'novel_denitrification_bubble.png')
out_pdf = os.path.join(WORK, 'novel_denitrification_bubble.pdf')
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print('wrote', out_png)
print('wrote', out_pdf)
