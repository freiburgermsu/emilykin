#!/usr/bin/env python3
"""
Gene abundance over time — one line per row of the gene_abundance heatmap.

Companion to 04_plot_heatmap.py. Collapses the per-MAG/per-sample RPKM matrix to
a community-level time series: for each KO functional row (the heatmap rows), the
abundance at each sample = sum of that row's per-sample RPKM across all MAGs
(RPKM is depth-normalised, so totals are comparable across samples / over time).

Outputs
  data/gene_abundance_over_time.tsv   spreadsheet: gene row × sample (+ summaries)
  gene_abundance_over_time.png / .pdf line figure: one line per heatmap row
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

# ── Heatmap rows, identical order/labels to 04_plot_heatmap.py ───────────────
# (stage, stage_color, [(row_id, short_label, full_label), ...])
ROW_GROUPS = [
    ('Nitrate reduction (NO₃→NO₂)', '#4C72B0',
     [('nar', 'nar', 'nar — membrane-bound nitrate reductase'),
      ('nap', 'nap', 'nap — periplasmic nitrate reductase')]),
    ('DNRA (NO₂→NH₄)', '#DD8452',
     [('nirB/D', 'nirB/D', 'nirB/D — NADH nitrite reductase (DNRA)')]),
    ('Denitrification (NO₂→NO)', '#55A868',
     [('nirS', 'nirS', 'nirS — cd₁ nitrite reductase'),
      ('nirK', 'nirK', 'nirK — Cu nitrite reductase')]),
    ('Denitrification (NO→N₂O)', '#C44E52',
     [('norC', 'norC', 'norC — nitric oxide reductase')]),
    ('Denitrification (N₂O→N₂)', '#8172B3',
     [('nosZ_cladeI', 'nosZ I', 'nosZ Clade I — N₂O reductase (typical)'),
      ('nosZ_cladeII', 'nosZ II', 'nosZ Clade II — N₂O reductase (atypical)')]),
]
# within-stage gene styling (1st gene solid/circle, 2nd dashed/square)
STYLES = [('-', 'o'), ('--', 's')]

# ── Aggregate: community-total RPKM per row per sample ───────────────────────
totals = {}        # row_id -> {sample -> summed RPKM}
nmags = {}         # row_id -> number of MAGs carrying the gene (RPKM>0 in any sample)
with open(os.path.join(DATA, 'gene_rpkm_per_sample_claded.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        g = r['ko_group']
        d = totals.setdefault(g, {s: 0.0 for s in SAMPLES})
        vals = {s: float(r[s]) for s in SAMPLES}
        for s in SAMPLES:
            d[s] += vals[s]
        if any(v > 0 for v in vals.values()):
            nmags[g] = nmags.get(g, 0) + 1

# flat, ordered row list mirroring the heatmap
rows = [(rid, short, full, stage, color)
        for stage, color, genes in ROW_GROUPS
        for rid, short, full in genes]

# ── Write the spreadsheet ───────────────────────────────────────────────────
tsv = os.path.join(DATA, 'gene_abundance_over_time.tsv')
with open(tsv, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko_group', 'gene', 'N_stage', *SAMPLES,
                'mean', 'max', 'n_MAGs_present'])
    for rid, short, full, stage, color in rows:
        vals = [totals[rid][s] for s in SAMPLES]
        w.writerow([rid, short, stage,
                    *[f'{v:.3f}' for v in vals],
                    f'{sum(vals)/len(vals):.3f}', f'{max(vals):.3f}',
                    nmags.get(rid, 0)])
print(f'Wrote {tsv}  ({len(rows)} gene rows × {len(SAMPLES)} samples)')

# ── Line figure: one line per heatmap row ───────────────────────────────────
x = list(range(1, len(SAMPLES) + 1))
fig, ax = plt.subplots(figsize=(8.4, 5.6), facecolor='white')

for stage, color, genes in ROW_GROUPS:
    for i, (rid, short, full) in enumerate(genes):
        ls, mk = STYLES[i] if i < len(STYLES) else ('-', 'o')
        y = [totals[rid][s] for s in SAMPLES]
        ax.plot(x, y, ls + mk, color=color, lw=1.9, ms=6, mew=0.8,
                mec='white', label=full, zorder=3)

ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(SAMPLES)
ax.set_xlim(0.85, len(SAMPLES) + 0.15)
ax.set_xlabel('Sample (CAN time series →)', fontsize=10)
ax.set_ylabel('Community gene abundance  (Σ RPKM over MAGs)', fontsize=10)
ax.set_title('Denitrification & DNRA gene abundance over time — CAN system',
             fontsize=12, fontweight='bold', pad=10)
ax.grid(True, which='major', axis='both', ls=':', lw=0.6, color='#cccccc', zorder=0)
ax.grid(True, which='minor', axis='y', ls=':', lw=0.4, color='#eeeeee', zorder=0)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)

# legend grouped by N-stage (titles act as group headers)
leg = ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5),
                fontsize=8, frameon=False, handlelength=2.6,
                title='N-transformation stage (heatmap rows)',
                title_fontsize=8.5, labelspacing=0.7)
leg._legend_box.align = 'left'

cap = ('Each line is one row of gene_abundance_figure (the heatmap). Abundance is the '
       'community total: per-sample RPKM summed across all 19 MAGs; RPKM is depth-'
       'normalised so values are comparable over time. Line colour = N-transformation '
       'stage; within a stage the second gene is dashed/square. Multi-gene rows '
       '(nar, nap, nirB/D, norC) use the per-MAG max RPKM over the row’s genes, as in '
       'the heatmap. Log y-axis. Source: data/gene_abundance_over_time.tsv.')
fig.text(0.012, -0.02, cap, ha='left', va='top', fontsize=7.0, family='serif',
         wrap=True)
fig.texts[-1]._get_wrap_line_width = lambda: 560

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(WORK, f'gene_abundance_over_time.{ext}'),
                dpi=200, bbox_inches='tight', facecolor='white')
print('Saved gene_abundance_over_time.png / .pdf')
plt.close()
