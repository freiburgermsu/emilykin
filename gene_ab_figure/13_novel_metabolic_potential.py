#!/usr/bin/env python3
"""
Metabolic-potential matrix of the NOVEL bacteria — Figure-S4-style companion to
the denitrification bubble plot (10_novel_denitrification_bubble.py).

Figure S4 of the reference (Wang/Gao/Wells, emi15486) shows, per dominant MAG,
the major carbon / electron-transport / EPS / P-metabolic pathways.  This is the
faithful *metagenomic* analog for this project's 7 species-novel MAGs: a
gene × MAG presence/copy-number matrix spanning the same functional axes that the
supinfo narrative uses to differentiate organism roles —

  carbon-source activation (acetate, propionate) · storage polymers (glycogen,
  PHA) · GAO gluconeogenesis (gapA/pckA/ppdk) · polyphosphate + phosphate (PAO) ·
  type IV pili (biofilm/granule aggregation) · vitamin-B12 uptake & biosynthesis ·
  biotin (cofactor) biosynthesis · secreted glycoside hydrolases (scavenging).

  rows    = marker genes, grouped/bracketed by functional category
  columns = the 7 novel MAGs, ordered as in the denitrification figure
            (most-complete denitrifier left), phylum header + N-cycle role strip
  cell    = gene copies per genome (discrete Blues, 1/2/3+); blank = absent

Data (./data): novel_functional_ko_long.tsv (08-style copy numbers from
11_novel_functional_annotation.py), gene_copy_per_mag.tsv, selected_mags.tsv,
taxonomy_labels.tsv, ../gtdbtk.bac120.summary.tsv.  matplotlib + numpy only.

Genomic potential, NOT expression (no metatranscriptome). CAZy/MEROPS/eggNOG deep
annotation was not run — glycoside-hydrolase rows are KO-level markers only and
secreted-protease/EPS inventories are out of scope.
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')
ROOT = os.path.dirname(WORK)

def load_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))

# ── functional copy numbers (category, KO, gene, MAG, copies) ────────────────
func = load_tsv(os.path.join(DATA, 'novel_functional_ko_long.tsv'))
copies = {}                                   # (gene, mag) -> copies
order_cat, order_gene = [], []                # preserve PANEL order
for r in func:
    copies[(r['gene'], r['MAG'])] = int(r['copies'])
    if r['category'] not in order_cat:
        order_cat.append(r['category'])
    if (r['category'], r['gene']) not in order_gene:
        order_gene.append((r['category'], r['gene']))

# ── column order: same as the denitrification bubble (most-complete left) ────
cp = {r['MAG']: r for r in load_tsv(os.path.join(DATA, 'gene_copy_per_mag.tsv'))}
DENIT_MEMBERS = ['narG', 'narH', 'narI', 'napA', 'napB', 'nirB', 'nirD', 'nirS',
                 'nirK', 'norB', 'norC', 'nosZ_cladeI', 'nosZ_cladeII']
iterid = {r['MAG']: r['iterativeID'] for r in cp.values()}
mags = sorted({r['MAG'] for r in func},
              key=lambda m: (-sum(int(float(cp[m].get(c, 0))) for c in DENIT_MEMBERS), iterid[m]))

# ── taxonomy + role for the column headers ───────────────────────────────────
role = {r['MAG']: r['classification'] for r in load_tsv(os.path.join(DATA, 'selected_mags.tsv'))}
phylum = {}
for r in load_tsv(os.path.join(ROOT, 'gtdbtk.bac120.summary.tsv')):
    d = dict(x.split('__', 1) for x in r['classification'].split(';'))
    p = d.get('p', '').strip()
    phylum[r['user_genome']] = {'Pseudomonadota': 'Proteobacteria'}.get(p, p)

CAT_COLORS = {
    'PAO': '#7B3294', 'Denitrifier/PAO': '#1B9E77', 'Denitrifier': '#E6AB02',
    'GAO': '#D95F02', 'GAO/PAO': '#E7298A', 'Denitrifier/GAO': '#A6761D',
    'Denitrifier/GAO/PAO': '#66A61E', 'Other': '#999999',
}

# ── matrix ───────────────────────────────────────────────────────────────────
genes = [g for _, g in order_gene]
nrow, ncol = len(genes), len(mags)
M = np.full((nrow, ncol), 0, dtype=int)
for i, (_, g) in enumerate(order_gene):
    for j, m in enumerate(mags):
        M[i, j] = copies.get((g, m), 0)

# category -> contiguous row span (PANEL groups genes by category already)
cat_span = {}
for i, (c, _) in enumerate(order_gene):
    cat_span.setdefault(c, [i, i])[1] = i

# ── plot ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(11.0, 14.2))

cmap = plt.get_cmap('Blues')
disc = {1: cmap(0.38), 2: cmap(0.62), 3: cmap(0.88)}

# faint background + filled cells
for i in range(nrow):
    for j in range(ncol):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor='#f5f5f5',
                               edgecolor='white', lw=0.8, zorder=0))
        cn = M[i, j]
        if cn > 0:
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                   facecolor=disc[min(cn, 3)], edgecolor='white',
                                   lw=0.8, zorder=1))
            txt = '3+' if cn >= 3 else str(cn)
            ax.text(j, i, txt, ha='center', va='center', fontsize=7.5, zorder=2,
                    color='white' if cn >= 2 else '#1f3b57')

ax.set_xlim(-0.5, ncol - 0.5)
ax.set_ylim(-0.5, nrow - 0.5)
ax.invert_yaxis()
ax.set_xticks(range(ncol))
ax.set_xticklabels([iterid[m] for m in mags], rotation=40, ha='right', fontsize=8.5)
ax.set_yticks(range(nrow))
ax.set_yticklabels([g for g in genes], fontsize=7.8, style='italic')
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)

# ── left category brackets + labels ──────────────────────────────────────────
xb = -2.05
for c in order_cat:
    r0, r1 = cat_span[c]
    ax.plot([xb, xb], [r0 - 0.45, r1 + 0.45], color='#444', lw=1.4,
            clip_on=False, zorder=5)
    for yy in (r0 - 0.45, r1 + 0.45):
        ax.plot([xb, xb + 0.10], [yy, yy], color='#444', lw=1.4, clip_on=False, zorder=5)
    ax.text(xb - 0.12, (r0 + r1) / 2, c, va='center', ha='right',
            fontsize=8.6, fontweight='bold', clip_on=False)

# ── top phylum headers + N-cycle role strip ──────────────────────────────────
strip_y = -0.5
for j, m in enumerate(mags):
    ax.add_patch(Rectangle((j - 0.46, strip_y - 0.95), 0.92, 0.55,
                           color=CAT_COLORS.get(role.get(m, 'Other'), '#999'),
                           clip_on=False, zorder=4))
    ax.text(j, strip_y - 1.08, phylum.get(m, ''), rotation=35, ha='left',
            va='bottom', fontsize=7.2, color='#555', clip_on=False)

# ── legends ──────────────────────────────────────────────────────────────────
copy_handles = [Line2D([], [], marker='s', linestyle='', markersize=11,
                       markerfacecolor=disc[c], markeredgecolor='white')
                for c in (1, 2, 3)]
leg1 = ax.legend(copy_handles, ['1', '2', '3+'], title='gene copies\nper genome',
                 loc='upper left', bbox_to_anchor=(1.02, 1.00), frameon=False,
                 fontsize=8.5, title_fontsize=8.8, labelspacing=0.6)
ax.add_artist(leg1)

present_roles = []
for m in mags:
    r = role.get(m, 'Other')
    if r not in present_roles:
        present_roles.append(r)
role_handles = [Rectangle((0, 0), 1, 1, color=CAT_COLORS.get(r, '#999')) for r in present_roles]
leg2 = ax.legend(role_handles, present_roles, title='N-cycle role', loc='upper left',
                 bbox_to_anchor=(1.02, 0.78), frameon=False, fontsize=8.3,
                 title_fontsize=8.8, handlelength=1.1, handleheight=1.1)
ax.add_artist(leg2)

# ── title + caption ──────────────────────────────────────────────────────────
fig.suptitle('Metabolic potential of the novel bacteria',
             x=0.02, y=0.992, ha='left', fontsize=14, fontweight='bold')
fig.text(0.02, 0.975,
         'species-novel MAGs (no GTDB species assignment) — metagenomic analog of Figure S4',
         ha='left', fontsize=9.3, color='#555')
cap = ('KOfam marker-gene inventory (pyhmmer, adaptive thresholds; Bakta proteomes) across the '
       'carbon / storage / phosphorus / pili / cofactor axes the reference supinfo uses to '
       'differentiate organism roles. Cell = gene copies per genome; blank = gene absent. Columns '
       'ordered as in the denitrification figure (most-complete denitrifier at left). Genomic '
       'potential, NOT expression. CAZy/MEROPS/eggNOG deep annotation was not run, so glycoside-'
       'hydrolase rows are KO-level markers only and secreted-protease / EPS inventories are out of scope.')
fig.text(0.02, 0.008, cap, ha='left', va='bottom', fontsize=7.0, color='#555', wrap=True)

fig.subplots_adjust(left=0.27, right=0.86, top=0.885, bottom=0.075)
out_png = os.path.join(WORK, 'novel_metabolic_potential.png')
out_pdf = os.path.join(WORK, 'novel_metabolic_potential.pdf')
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print('genes:', nrow, ' mags:', ncol)
print('column order:', [iterid[m] for m in mags])
print('wrote', out_png)
print('wrote', out_pdf)
