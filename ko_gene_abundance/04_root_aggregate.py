#!/usr/bin/env python3
"""
Condensed, organism-level KO abundance: sum the per-MAG KO RPKM across MAGs that
share a mag_iterativeID root (e.g. Ca_Accumulibacter.8 + Ca_Accumulibacter.23 +
... -> one "Ca_Accumulibacter" row). Because the 276-MAG community mapping
partitions a read to its single best MAG, same-organism bins (co-assembly +
per-sample) split the signal; summing by root recovers the organism-level total.

Root = mag_iterativeID with the trailing ".<n>" removed.

Inputs : ko_rpkm_combined.tsv, ko_rpkm_per_sample.tsv  (from 03)
Outputs: ko_rpkm_root_combined.tsv     (root x 42 KOs, combined RPKM)
         ko_rpkm_root_per_sample.tsv   (ko, root, CAN_1..CAN_5)
"""
import csv, os
from collections import defaultdict

WORK = os.path.dirname(os.path.abspath(__file__))
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

def root_of(iterid, mag):
    return iterid.rsplit('.', 1)[0] if iterid else mag

# ── combined per-MAG -> root ────────────────────────────────────────────────
rows = list(csv.DictReader(open(os.path.join(WORK, 'ko_rpkm_combined.tsv')), delimiter='\t'))
kos = [c for c in rows[0] if c.startswith('K') and c[1:].isdigit()]
root_comb = defaultdict(lambda: defaultdict(float))
root_mags = defaultdict(list)
root_genus = {}
for r in rows:
    root = root_of(r['mag_iterativeID'], r['MAG'])
    root_mags[root].append(r['MAG'])
    root_genus.setdefault(root, r['genus'])
    for ko in kos:
        root_comb[root][ko] += float(r[ko])

# order roots by total combined RPKM (most active organisms first)
roots = sorted(root_comb, key=lambda x: -sum(root_comb[x].values()))

with open(os.path.join(WORK, 'ko_rpkm_root_combined.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['root_iterativeID', 'n_MAGs', 'genus'] + kos)
    for root in roots:
        w.writerow([root, len(root_mags[root]), root_genus.get(root, '')] +
                   [f'{root_comb[root][ko]:.4f}' for ko in kos])

# ── per-sample per-MAG -> root ──────────────────────────────────────────────
root_ps = defaultdict(lambda: defaultdict(lambda: {s: 0.0 for s in SAMPLES}))
for r in csv.DictReader(open(os.path.join(WORK, 'ko_rpkm_per_sample.tsv')), delimiter='\t'):
    root = root_of(r['iterativeID'], r['mag'])
    for s in SAMPLES:
        root_ps[r['ko']][root][s] += float(r[s])

with open(os.path.join(WORK, 'ko_rpkm_root_per_sample.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko', 'root_iterativeID', 'n_MAGs', 'genus'] + SAMPLES)
    for ko in kos:
        for root in roots:
            v = root_ps[ko][root]
            if sum(v.values()) == 0 and root_comb[root][ko] == 0:
                continue
            w.writerow([ko, root, len(root_mags[root]), root_genus.get(root, '')] +
                       [f'{v[s]:.4f}' for s in SAMPLES])

print(f'{len(rows)} MAGs -> {len(roots)} organism roots')
print(f'roots with >1 MAG (signal-pooling): {sum(1 for r in roots if len(root_mags[r])>1)}')
print('wrote ko_rpkm_root_combined.tsv and ko_rpkm_root_per_sample.tsv')

# ── readable preview: denitrification KOs for the top denitrifier roots ──────
DENIT = [('K00370','narG'),('K00371','narH'),('K00374','narI'),('K02567','napA'),
         ('K02568','napB'),('K15864','nirS'),('K00368','nirK'),('K04561','norB'),
         ('K02305','norC'),('K00376','nosZ')]
denit_total = {root: sum(root_comb[root][k] for k, _ in DENIT) for root in roots}
top = sorted(roots, key=lambda r: -denit_total[r])[:18]
print('\nTop denitrifier organisms (combined RPKM; root = iterativeID base):')
print(f"  {'root':24s} {'nMAG':>4s} " + ' '.join(f'{n:>6s}' for _, n in DENIT))
for root in top:
    if denit_total[root] <= 0:
        continue
    print(f"  {root:24s} {len(root_mags[root]):>4d} " +
          ' '.join(f'{root_comb[root][k]:6.1f}' for k, _ in DENIT))
