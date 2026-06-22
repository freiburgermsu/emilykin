#!/usr/bin/env python3
"""
Aggregate per-sample read counts into a KO x MAG gene-abundance matrix (RPKM),
recreating the KO block of mag_abundance_summary.tsv as quantitative abundance
(vs the existing 0/1 presence).

per-gene RPKM(sample) = count / (gene_len_kb) / (mapped_reads_million)
KO abundance(MAG, sample) = SUM over the MAG's copies of that KO of their RPKM
  (total read abundance of the KO across its gene copies; with --secondary=no each
   read maps to one copy, so summing copies = total reads attributable to the KO).

Covers all 276 MAGs x 42 KOs from mag_abundance_summary.tsv (0 where the MAG has no
such gene). MAGs labelled by mag_iterativeID.

Outputs (ko_gene_abundance/):
  ko_rpkm_per_sample.tsv     long: ko, mag, iterativeID, genus, CAN_1..CAN_5
  ko_rpkm_combined.tsv       wide: MAG x KO (combined-across-samples RPKM)
"""
import csv, json, os
from collections import defaultdict

WORK = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(WORK)
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

# ── canonical MAG list, KO list, labels from mag_abundance_summary.tsv ───────
mags, kos, genus, iterid = [], [], {}, {}
with open(os.path.join(REPO, 'mag_abundance_summary.tsv')) as f:
    rdr = csv.DictReader(f, delimiter='\t')
    kos = [c for c in rdr.fieldnames if c.startswith('K') and c[1:].isdigit()]
    for r in rdr:
        m = r['MAG']; mags.append(m)
        genus[m] = r.get('Genus', ''); iterid[m] = r.get('mag_iterativeID', '')
print(f'{len(mags)} MAGs x {len(kos)} KOs (from mag_abundance_summary.tsv)')

# ── KO-gene loci: gene_id -> (len, mag, {kos}) ──────────────────────────────
glen, gmag, gkos = {}, {}, defaultdict(set)
with open(os.path.join(WORK, 'ko_genes.bed')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        gid = r['gene_id']
        glen[gid] = max(int(r['end']) - int(r['start']), 1)
        gmag[gid] = r['mag']
        gkos[gid].add(r['ko'])

# ── per-sample counts ───────────────────────────────────────────────────────
tot, count = {}, defaultdict(dict)
for s in SAMPLES:
    p = os.path.join(WORK, f'counts_{s}.json')
    if not os.path.exists(p):
        raise SystemExit(f'missing {p} — run 02_align_ko.sh first')
    d = json.load(open(p)); tot[s] = d['total_mapped']
    for gid, c in d['gene_counts'].items():
        count[gid][s] = c
    print(f'  {s}: {tot[s]:,} mapped reads')
tot_comb = sum(tot.values())

def rpkm(c, L, t):
    return (c / (L / 1000.0)) / (t / 1e6) if t and c else 0.0

# ── aggregate: (mag, ko) -> per-sample + combined RPKM (sum over copies) ─────
per = defaultdict(lambda: {s: 0.0 for s in SAMPLES})
comb = defaultdict(float)
for gid, kset in gkos.items():
    m, L = gmag[gid], glen[gid]
    cs = {s: count[gid].get(s, 0) for s in SAMPLES}
    tc = sum(cs.values())
    for ko in kset:
        for s in SAMPLES:
            per[(ko, m)][s] += rpkm(cs[s], L, tot[s])
        comb[(ko, m)] += rpkm(tc, L, tot_comb)

# ── write long per-sample ───────────────────────────────────────────────────
with open(os.path.join(WORK, 'ko_rpkm_per_sample.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko', 'mag', 'iterativeID', 'genus'] + SAMPLES)
    for ko in kos:
        for m in mags:
            v = per.get((ko, m))
            if v is None and comb.get((ko, m), 0) == 0:
                row_vals = [0.0] * len(SAMPLES)
            else:
                row_vals = [v[s] if v else 0.0 for s in SAMPLES]
            w.writerow([ko, m, iterid.get(m, ''), genus.get(m, '')] +
                       [f'{x:.4f}' for x in row_vals])

# ── write wide combined MAG x KO (mirrors mag_abundance_summary KO block) ────
with open(os.path.join(WORK, 'ko_rpkm_combined.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['MAG', 'mag_iterativeID', 'genus'] + kos)
    for m in mags:
        w.writerow([m, iterid.get(m, ''), genus.get(m, '')] +
                   [f'{comb.get((ko, m), 0.0):.4f}' for ko in kos])

# ── summary ─────────────────────────────────────────────────────────────────
present = sum(1 for k in comb if comb[k] > 0)
print(f'\nwrote ko_rpkm_per_sample.tsv and ko_rpkm_combined.tsv')
print(f'{present} (KO,MAG) cells with abundance > 0')
print('per-KO: MAGs with signal / max combined RPKM:')
for ko in kos:
    vals = [comb.get((ko, m), 0.0) for m in mags]
    nz = sum(1 for v in vals if v > 0)
    print(f'  {ko}: {nz:3d} MAGs, max={max(vals):.1f}')
