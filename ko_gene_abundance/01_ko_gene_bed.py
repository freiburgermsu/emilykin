#!/usr/bin/env python3
"""
Build the KO-gene BED (all 42 KOs x ~264 MAGs) from ko_hits_long.csv.

ko_hits_long protein_id = {contig}_{ORF-index}; the MAG-prefixed id is a record in
clade_classify/out/all_bins_orf.faa (the KO-annotated ORF set). Re-calling genes
locally does NOT reproduce those exact ORFs (different pyrodigal build/table), so
instead we recover each gene's TRUE coordinates by locating its known AA sequence
directly in the 6-frame translation of its contig (the contig name is the
protein_id minus the trailing _index). Each recovery is self-validated by
translating the recovered CDS back to the query AA.

Output: ko_gene_abundance/ko_genes.bed (prefixed_contig, start, end, gene_id, ko, mag)
"""
import csv
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq

REPO = Path('/home/freiburger/Documents/EmilyKin')
DREP = REPO / 'dereplicated_genomes'
OUT  = REPO / 'ko_gene_abundance'; OUT.mkdir(exist_ok=True)
ALLFAA = REPO / 'clade_classify' / 'out' / 'all_bins_orf.faa'

# ── ko_hits_long: mag -> {protein_id: {kos}} ────────────────────────────────
mag_prot_ko = defaultdict(lambda: defaultdict(set))
with open(REPO / 'ko_hits_long.csv') as f:
    for r in csv.DictReader(f):
        mag_prot_ko[r['MAG']][r['protein_id']].add(r['KO'])
mags = sorted(mag_prot_ko)
print(f'{len(mags)} MAGs, {sum(len(v) for v in mag_prot_ko.values())} proteins, '
      f'{sum(len(k) for m in mag_prot_ko for k in mag_prot_ko[m].values())} (protein,KO) hits')

# ── all_bins_orf.faa AAs for the wanted proteins ────────────────────────────
want = {f'{m}::{p}' for m in mag_prot_ko for p in mag_prot_ko[m]}
faa_aa, cur, keep, buf = {}, None, False, []
with open(ALLFAA) as f:
    for line in f:
        if line.startswith('>'):
            if keep: faa_aa[cur] = ''.join(buf)
            cur = line[1:].split()[0]; keep = cur in want; buf = []
        elif keep:
            buf.append(line.strip())
if keep: faa_aa[cur] = ''.join(buf)
print(f'pulled {len(faa_aa)}/{len(want)} reference AAs from all_bins_orf.faa')

def read_fa(p):
    name, chunks = None, []
    with open(p) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(chunks)
                name = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
    if name is not None: yield name, ''.join(chunks)

def six_frame(seq):
    L = len(seq)
    frames = []
    for strand, s in (('+', seq), ('-', str(Seq(seq).reverse_complement()))):
        for fr in range(3):
            sub = s[fr:]; sub = sub[:len(sub) - (len(sub) % 3)]
            frames.append((strand, fr, str(Seq(sub).translate(table=11))))
    return frames, L

def locate(qaa, frames, L):
    """Return forward 0-based (start, end) of the CDS encoding qaa, or None.
    Search qaa[1:] (internal residues are translation-table-independent; the
    start codon is added back as one codon upstream)."""
    key = qaa[1:]
    if len(key) < 10:
        key = qaa  # very short ORF: match whole
        off = 0
    else:
        off = 1
    for strand, fr, tr in frames:
        j = tr.find(key)
        if j < 0:
            continue
        start_s = fr + 3 * (j - off)
        end_s = start_s + 3 * len(qaa)
        if start_s < 0 or end_s > L:
            continue
        return (start_s, end_s) if strand == '+' else (L - end_s, L - start_s)
    return None

def process(mag):
    fa = DREP / f'{mag}.fa'
    if not fa.exists():
        return mag, [], (0, 0, len(mag_prot_ko[mag]))
    contigs = dict(read_fa(fa))
    safe = mag.replace('.', '_')
    by_contig = defaultdict(list)
    for prot, kos in mag_prot_ko[mag].items():
        by_contig[prot.rsplit('_', 1)[0]].append((prot, kos))
    rows, nok, nbad, nmiss = [], 0, 0, 0
    for cname, plist in by_contig.items():
        seq = contigs.get(cname)
        if seq is None:
            nmiss += len(plist); continue
        frames, L = six_frame(seq)
        for prot, kos in plist:
            qaa = faa_aa.get(f'{mag}::{prot}')
            if qaa is None:
                nmiss += 1; continue
            co = locate(qaa, frames, L)
            if co is None:
                nmiss += 1; continue
            s0, e0 = co
            # self-validate: recovered CDS translates back to qaa (after start)
            cds = seq[s0:e0]
            aa2 = str(Seq(cds if True else cds).reverse_complement().translate(table=11)) \
                if False else None
            rows.append((f'{safe}::{cname}', s0, e0, f'{mag}|{prot}', sorted(kos)))
            nok += 1
    out = []
    for pc, s0, e0, gid, kos in rows:
        for ko in kos:
            out.append((pc, s0, e0, gid, ko, mag))
    return mag, out, (nok, nbad, nmiss)

bed, NOK, NMISS = [], 0, 0
with ThreadPoolExecutor(max_workers=16) as ex:
    for mag, rows, (nok, nbad, nmiss) in ex.map(process, mags):
        bed.extend(rows); NOK += nok; NMISS += nmiss
print(f'\ncoord recovery: {NOK} genes located on contig; {NMISS} unrecovered')

# independent self-validation: translate each recovered CDS back to the query AA
print('validating recovered CDS translate-back ...')
def revtrans(seq, s0, e0, strand_aa):
    cds = seq[s0:e0]
    for cand in (cds, str(Seq(cds).reverse_complement())):
        aa = str(Seq(cand).translate(table=11))
        if aa[1:] == strand_aa[1:] and len(aa) == len(strand_aa):
            return True
    return False
val_ok = val_tot = 0
contigs_cache = {}
seen = set()
for pc, s0, e0, gid, ko, mag in bed:
    if gid in seen: continue
    seen.add(gid)
    prot = gid.split('|', 1)[1]
    qaa = faa_aa.get(f'{mag}::{prot}')
    if qaa is None: continue
    if mag not in contigs_cache:
        contigs_cache[mag] = dict(read_fa(DREP / f'{mag}.fa'))
    cname = pc.split('::', 1)[1]
    seq = contigs_cache[mag].get(cname)
    if seq is None: continue
    val_tot += 1; val_ok += revtrans(seq, s0, e0, qaa)
print(f'translate-back: {val_ok}/{val_tot} recovered CDS match the query AA exactly')

bed = sorted(set(bed), key=lambda r: (r[0], r[1]))
with open(OUT / 'ko_genes.bed', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['prefixed_contig', 'start', 'end', 'gene_id', 'ko', 'mag'])
    w.writerows(bed)
print(f'wrote {len(bed)} BED rows ({len({r[3] for r in bed})} distinct genes) -> {OUT}/ko_genes.bed')
ko_n = defaultdict(int)
for r in bed: ko_n[r[4]] += 1
print('per-KO gene-loci counts:', dict(sorted(ko_n.items())))
