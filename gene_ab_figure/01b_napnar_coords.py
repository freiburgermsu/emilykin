#!/usr/bin/env python3
"""
Add nar/nap nitrate-reductase genes to the gene-abundance target set.

The broad 42-KO KOfam search (kofamscan/hmmsearch.tblout) already scored nar/nap
in every MAG; this step (a) selects the passing nar/nap proteins in the figure's
selected MAGs, and (b) recovers their contig coordinates. bakta GFF3 (the original
coord source in 01_select_mags_and_genes.py) is not present locally, so coords are
recovered by re-calling each MAG with pyrodigal (single mode = what bakta uses) and
matching the bakta protein to its ORF. The method is VALIDATED here by re-deriving
the coordinates of existing target_genes.bed entries and checking they agree.

KO subtypes (max-RPKM policy applied later, like nirB/D and norC):
  nar : K00370 narG, K00371 narH, K00374 narI   (membrane-bound; NO3->NO2)
  nap : K02567 napA, K02568 napB                 (periplasmic;    NO3->NO2)

Outputs:
  data/napnar_genes.bed      BED rows to append to target_genes.bed
  (prints a validation table)
"""
import csv, sys
from pathlib import Path
import pyrodigal

WORK = Path('/home/freiburger/Documents/EmilyKin/gene_ab_figure')
DATA = WORK / 'data'
REPO = Path('/home/freiburger/Documents/EmilyKin')
KOFAM = REPO / 'kofamscan'
DREP = REPO / 'dereplicated_genomes'

KO = {'K00370':'narG','K00371':'narH','K00374':'narI','K02567':'napA','K02568':'napB'}

# ── thresholds + selected MAGs ──────────────────────────────────────────────
thr = {}
with open(REPO/'meta'/'mag_gene_ab'/'kofam'/'target_ko_list.txt') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['knum'] in KO:
            thr[r['knum']] = (float(r['threshold']), r['score_type'])
selected = [l.strip() for l in open(DATA/'selected_mag_list.txt') if l.strip()]
selset = set(selected)

# ── nar/nap hits passing threshold, in selected MAGs ────────────────────────
napnar = []   # (mag, prot, ko, gene)
with open(KOFAM/'hmmsearch.tblout') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.split()
        if len(p) < 10: continue
        target, ko = p[0], p[2]
        if ko not in KO or '|' not in target: continue
        mag, prot = target.split('|', 1)
        if mag not in selset: continue
        score = float(p[8]) if thr[ko][1] == 'domain' else float(p[5])
        if score >= thr[ko][0]:
            napnar.append((mag, prot, ko, KO[ko]))
nap_mags = sorted({m for m, *_ in napnar})
print(f'nar/nap genes to place: {len(napnar)} across {len(nap_mags)} MAGs')

# ── existing BED genes in the same MAGs → validation set + their coords ─────
bed_coords = {}   # (mag, locus) -> (prefixed_contig, start0, end)
with open(DATA/'target_genes.bed') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        m = r['mag']
        if m in set(nap_mags):
            locus = r['gene_id'].split('|', 1)[1]
            bed_coords[(m, locus)] = (r['prefixed_contig'], int(r['start']), int(r['end']))

want_prots = {f'{m}|{prot}' for m, prot, _, _ in napnar} | {f'{m}|{loc}' for (m, loc) in bed_coords}

# ── pull AA sequences for wanted proteins from all_mags.faa ──────────────────
aa = {}
cur, keep, buf = None, False, []
with open(KOFAM/'all_mags.faa') as f:
    for line in f:
        if line.startswith('>'):
            if keep: aa[cur] = ''.join(buf)
            cur = line[1:].split()[0]; keep = cur in want_prots; buf = []
        elif keep:
            buf.append(line.strip())
if keep: aa[cur] = ''.join(buf)
print(f'pulled {len(aa)}/{len(want_prots)} protein sequences from all_mags.faa')

# ── per-MAG pyrodigal (single mode) → AA-indexed coordinate table ────────────
def read_fa(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(chunks)
                name = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
    if name is not None: yield name, ''.join(chunks)

def orf_index(mag):
    """Return (exact: aa->coords, ctail: last30aa->[coords]) for one MAG."""
    contigs = list(read_fa(DREP/f'{mag}.fa'))
    gf = pyrodigal.GeneFinder()
    try:
        gf.train(*[s.encode() for _, s in contigs])
    except Exception:
        gf = pyrodigal.GeneFinder(meta=True)
    exact, ctail = {}, {}
    safe = mag.replace('.', '_')
    for cname, seq in contigs:
        for g in gf.find_genes(seq.encode()):
            a = g.translate().rstrip('*')
            coords = (f'{safe}::{cname}', g.begin - 1, g.end, '+' if g.strand == 1 else '-')
            exact.setdefault(a, coords)
            ctail.setdefault(a[-30:], []).append((a, coords))
    return exact, ctail

def match(prot_aa, exact, ctail):
    if prot_aa in exact:
        return exact[prot_aa], 'exact'
    # tolerant: identical C-terminus (callers agree on stop), len within 15 aa
    for a, coords in ctail.get(prot_aa[-30:], []):
        if abs(len(a) - len(prot_aa)) <= 15:
            return coords, f'ctail(Δlen={len(a)-len(prot_aa)})'
    return None, 'NO_MATCH'

# ── validate on existing BED genes, then place nar/nap ───────────────────────
print('\n=== VALIDATION: re-derive existing BED gene coords via pyrodigal ===')
nok = nbad = 0
napnar_rows = []
for mag in nap_mags:
    exact, ctail = orf_index(mag)
    # validation
    for (m, loc), (pc, s0, e0) in bed_coords.items():
        if m != mag: continue
        pid = f'{m}|{loc}'
        if pid not in aa: continue
        coords, how = match(aa[pid], exact, ctail)
        if coords is None:
            print(f'  [val] {pid:32s} NO pyrodigal match'); nbad += 1; continue
        pc2, s2, e2 = coords[0], coords[1], coords[2]
        ok = (pc2 == pc) and (abs(s2 - s0) <= 30) and (abs(e2 - e0) <= 30)
        nok += ok; nbad += (not ok)
        if not ok:
            print(f'  [val] {pid:32s} MISMATCH bed=({pc},{s0},{e0}) pyro=({pc2},{s2},{e2}) [{how}]')
    # place nar/nap
    for (m, prot, ko, gene) in napnar:
        if m != mag: continue
        pid = f'{m}|{prot}'
        coords, how = match(aa.get(pid, ''), exact, ctail)
        if coords is None:
            print(f'  [napnar] {pid} ({gene}) NO MATCH — needs manual coords'); continue
        pc, s0, e0, strand = coords
        napnar_rows.append((pc, s0, e0, f'{m}|{prot}', ko, m, gene, how))
print(f'validation: {nok} agree, {nbad} disagree (≤30 bp tolerance, same contig)')

# ── write BED additions ─────────────────────────────────────────────────────
out = DATA / 'napnar_genes.bed'
with open(out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    for pc, s0, e0, gid, ko, mag, gene, how in sorted(napnar_rows):
        w.writerow([pc, s0, e0, gid, ko, mag])
print(f'\nwrote {len(napnar_rows)} nar/nap BED rows -> {out}')
print('\n=== nar/nap genes placed ===')
for pc, s0, e0, gid, ko, mag, gene, how in sorted(napnar_rows):
    print(f'  {gene:5s} {ko}  {gid:30s} {pc}:{s0}-{e0}  ({e0-s0} bp) [{how}]')
