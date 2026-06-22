#!/usr/bin/env python3
"""
Aggregate per-sample counts_<sample>.json into the clade-split RPKM matrix the
heatmap consumes (data/gene_rpkm_per_sample_claded.tsv), now including the
nitrate-reductase rows.

Self-contained replacement for the (stale) 03_aggregate + 05 §6 chain: reads
total_mapped + per-gene counts directly from counts_<sample>.json, so it needs no
total-read recovery. RPKM = count/(len_kb)/(mapped_million). Each figure row's
cell per MAG/sample = MAX RPKM over the row's KO subtypes (same policy as nirB/D
and norC):

  nar          K00370 narG / K00371 narH / K00374 narI   (NO3->NO2, membrane)
  nap          K02567 napA / K02568 napB                 (NO3->NO2, periplasmic)
  nirB/D       K00362 / K00363
  nirS         K15864
  nirK         K00368
  norC         K02305 norC / K04561 norB
  nosZ_cladeI  K00376 where nosz_clades.tsv == I
  nosZ_cladeII K00376 where nosz_clades.tsv == II

Validates that the six pre-existing rows reproduce the committed claded.tsv.
"""
import csv, json, os
from collections import defaultdict

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

# row -> set of KOs (nosZ handled separately via clade)
ROW_KOS = {
    'nar':    {'K00370', 'K00371', 'K00374'},
    'nap':    {'K02567', 'K02568'},
    'nirB/D': {'K00362', 'K00363'},
    'nirS':   {'K15864'},
    'nirK':   {'K00368'},
    'norC':   {'K02305', 'K04561'},
}
OUT_ROWS = ['nar', 'nap', 'nirB/D', 'nirS', 'nirK', 'norC', 'nosZ_cladeI', 'nosZ_cladeII']

# ── target genes: gene_id -> (len, mag, {kos}) ──────────────────────────────
glen, gmag, gkos = {}, {}, defaultdict(set)
mags_in_bed = []
with open(os.path.join(DATA, 'target_genes.bed')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        gid = r['gene_id']
        glen[gid] = max(int(r['end']) - int(r['start']), 1)
        gmag[gid] = r['mag']
        gkos[gid].add(r['ko'])
        if r['mag'] not in mags_in_bed:
            mags_in_bed.append(r['mag'])
mag_list = sorted(mags_in_bed)

# ── nosZ clade per gene ─────────────────────────────────────────────────────
clade_of = {}
with open(os.path.join(DATA, 'nosz_clades.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        clade_of[r['gene_id']] = r['clade']

# ── per-sample counts ───────────────────────────────────────────────────────
tot = {}
count = defaultdict(dict)   # gid -> sample -> count
for s in SAMPLES:
    p = os.path.join(DATA, f'counts_{s}.json')
    if not os.path.exists(p):
        raise SystemExit(f'missing {p} — run 02_align_local.sh first')
    d = json.load(open(p))
    tot[s] = d['total_mapped']
    for gid, c in d['gene_counts'].items():
        count[gid][s] = c
    print(f'  {s}: {tot[s]:,} mapped reads')

def rpkm(gid, s):
    c = count[gid].get(s, 0)
    return (c / (glen[gid] / 1000.0)) / (tot[s] / 1e6) if tot[s] and c else 0.0

def row_match(rowname, gid):
    kos = gkos[gid]
    if rowname == 'nosZ_cladeI':
        return 'K00376' in kos and clade_of.get(gid) == 'I'
    if rowname == 'nosZ_cladeII':
        return 'K00376' in kos and clade_of.get(gid) == 'II'
    return bool(kos & ROW_KOS[rowname])

# ── build matrix (max over matching genes per MAG/sample) ───────────────────
out = [['ko_group', 'mag'] + SAMPLES]
for rowname in OUT_ROWS:
    for mag in mag_list:
        gids = [g for g in glen if gmag[g] == mag and row_match(rowname, g)]
        vals = [max((rpkm(g, s) for g in gids), default=0.0) for s in SAMPLES]
        out.append([rowname, mag] + [f'{v:.4f}' for v in vals])

# ── validate the 6 pre-existing rows against committed claded.tsv ───────────
prev = {}
pp = os.path.join(DATA, 'gene_rpkm_per_sample_claded.tsv')
if os.path.exists(pp):
    for r in csv.DictReader(open(pp), delimiter='\t'):
        prev[(r['ko_group'], r['mag'])] = {s: float(r[s]) for s in SAMPLES}
maxdiff = 0.0; ndiff = 0
for row in out[1:]:
    key = (row[0], row[1])
    if key in prev:
        for i, s in enumerate(SAMPLES):
            d = abs(float(row[2 + i]) - prev[key][s])
            if d > 0.01:
                ndiff += 1
            maxdiff = max(maxdiff, d)
print(f'\nvalidation vs committed claded.tsv (existing rows): max abs diff = {maxdiff:.4f}, cells>0.01: {ndiff}')

with open(pp, 'w', newline='') as f:
    csv.writer(f, delimiter='\t').writerows(out)
print(f'wrote {pp}  ({len(OUT_ROWS)} rows x {len(mag_list)} MAGs)')

# ── provenance: read_counts_raw.tsv ─────────────────────────────────────────
with open(os.path.join(DATA, 'read_counts_raw.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['gene_id', 'ko', 'mag', 'gene_len_bp'] + SAMPLES + ['combined_count'])
    for gid in sorted(glen):
        cs = [count[gid].get(s, 0) for s in SAMPLES]
        w.writerow([gid, '/'.join(sorted(gkos[gid])), gmag[gid], glen[gid]] + cs + [sum(cs)])

# ── per-row presence summary ────────────────────────────────────────────────
print('\nMAGs with signal per row:')
for rowname in OUT_ROWS:
    n = sum(1 for row in out[1:] if row[0] == rowname and any(float(x) > 0 for x in row[2:]))
    print(f'  {rowname:13s}: {n}')
