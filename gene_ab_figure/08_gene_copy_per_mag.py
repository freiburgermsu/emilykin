#!/usr/bin/env python3
"""
Genomic gene-copy inventory: number of copies of each denitrification/DNRA gene
(by KO) in each selected MAG, with nosZ (K00376) separated into Clade I and
Clade II. Reference-independent — counted directly from the bakta-derived gene
loci in data/target_genes.bed (one row per detected KO locus), so it is unaffected
by read-mapping reference choice.

Outputs:
  data/gene_copy_per_mag.tsv  (MAG x gene copy-count matrix)
  data/gene_copy_per_mag.md   (same, GitHub-flavoured markdown)
"""
import csv, os
from collections import defaultdict

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')

# gene columns in pathway order (NO3 -> NO2 -> NO -> N2O -> N2); nosZ split by clade
GENE_KOS = [('narG', 'K00370'), ('narH', 'K00371'), ('narI', 'K00374'),
            ('napA', 'K02567'), ('napB', 'K02568'),
            ('nirB', 'K00362'), ('nirD', 'K00363'),
            ('nirS', 'K15864'), ('nirK', 'K00368'),
            ('norB', 'K04561'), ('norC', 'K02305')]
NOSZ_KO = 'K00376'

# ── copies: distinct gene loci per (mag, ko) from the BED ────────────────────
copies = defaultdict(lambda: defaultdict(set))   # mag -> ko -> {gene_id}
nosz_genes = defaultdict(list)                    # mag -> [gene_id]  (K00376)
with open(os.path.join(DATA, 'target_genes.bed')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        copies[r['mag']][r['ko']].add(r['gene_id'])
        if r['ko'] == NOSZ_KO:
            nosz_genes[r['mag']].append(r['gene_id'])

# ── nosZ clade + Clade-II sub-clade per gene ────────────────────────────────
clade_of = {}
with open(os.path.join(DATA, 'nosz_clades.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        clade_of[r['gene_id']] = r['clade']
subclade_of = {}
scp = os.path.join(DATA, 'nosz_cladeII_subclades.tsv')
if os.path.exists(scp):
    with open(scp) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            subclade_of[r['gene_id']] = r['subclade']

# ── MAG labels ──────────────────────────────────────────────────────────────
genus, iterid = {}, {}
with open(os.path.join(DATA, 'taxonomy_labels.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        genus[r['MAG']] = r.get('Genus', '')
        iterid[r['MAG']] = r.get('mag_iterativeID', '')
mags = [l.strip() for l in open(os.path.join(DATA, 'selected_mag_list.txt')) if l.strip()]
mags.sort()

# ── ppk1 classification (Accumulibacter MAGs only; ppk1_classify/) ───────────
ppk1_type, ppk1_species = {}, {}
ppk1_path = os.path.join(WORK, '..', 'ppk1_classify', 'ppk1_classification.tsv')
if os.path.exists(ppk1_path):
    with open(ppk1_path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            ppk1_type[r['MAG']] = r.get('ppk1_type', '')
            ppk1_species[r['MAG']] = r.get('proposed_species', '')

# ── build rows ──────────────────────────────────────────────────────────────
COLS = [g for g, _ in GENE_KOS] + ['nosZ_cladeI', 'nosZ_cladeII']
rows = []
for m in mags:
    cnt = {g: len(copies[m].get(ko, ())) for g, ko in GENE_KOS}
    cI = sum(1 for gid in nosz_genes[m] if clade_of.get(gid) == 'I')
    cII = sum(1 for gid in nosz_genes[m] if clade_of.get(gid) == 'II')
    cnt['nosZ_cladeI'], cnt['nosZ_cladeII'] = cI, cII
    sub = ';'.join(sorted({subclade_of[g] for g in nosz_genes[m]
                           if clade_of.get(g) == 'II' and g in subclade_of})) or '-'
    total = sum(cnt[c] for c in COLS)
    rows.append({'MAG': m, 'iterativeID': iterid.get(m, ''), 'genus': genus.get(m, ''),
                 **cnt, 'total': total, 'nosZ_II_subclade': sub,
                 'ppk1_type': ppk1_type.get(m, ''), 'ppk1_species': ppk1_species.get(m, '')})

# ── write TSV ────────────────────────────────────────────────────────────────
fields = (['MAG', 'iterativeID', 'genus'] + COLS
          + ['total', 'nosZ_II_subclade', 'ppk1_type', 'ppk1_species'])
with open(os.path.join(DATA, 'gene_copy_per_mag.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
    tot = {'MAG': 'TOTAL', 'iterativeID': '', 'genus': '',
           **{c: sum(r[c] for r in rows) for c in COLS},
           'total': sum(r['total'] for r in rows), 'nosZ_II_subclade': '',
           'ppk1_type': '', 'ppk1_species': ''}
    w.writerow(tot)

# ── write markdown ───────────────────────────────────────────────────────────
def md_table(fields, rows, tot):
    hdr = '| ' + ' | '.join(fields) + ' |'
    sep = '|' + '|'.join('---' if i < 3 else ':-:' for i in range(len(fields))) + '|'
    out = [hdr, sep]
    for r in rows:
        out.append('| ' + ' | '.join(str(r[c]) if r[c] != 0 else '·' for c in fields) + ' |')
    out.append('| ' + ' | '.join(str(tot.get(c, '')) if tot.get(c, '') != 0 else '·' for c in fields) + ' |')
    return '\n'.join(out)

legend = ('Gene→KO: narG K00370, narH K00371, narI K00374, napA K02567, napB K02568, '
          'nirB K00362, nirD K00363, nirS K15864, nirK K00368, norB K04561, norC K02305, '
          'nosZ K00376 (split by Clade I/II by tree placement, reconciled — '
          'clade_classify/partD; subclade via nhmmer vs NosZREF). '
          'ppk1_type/ppk1_species: Accumulibacter ppk1 clade + proposed species '
          '(Petriglieri 2022; ppk1_classify/). '
          'Copies = distinct gene loci per MAG (genomic inventory; · = 0).')
with open(os.path.join(DATA, 'gene_copy_per_mag.md'), 'w') as f:
    f.write('# Denitrification / DNRA gene copies per MAG (genomic inventory)\n\n')
    f.write(legend + '\n\n')
    f.write(md_table(fields, rows, tot) + '\n')

# ── console ──────────────────────────────────────────────────────────────────
print(md_table(fields, rows, tot))
print('\n' + legend)
print(f'\nwrote {DATA}/gene_copy_per_mag.tsv and .md')
