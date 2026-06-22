#!/usr/bin/env python3
"""
Assign Clade II sub-clade (NosZREF 1577 A-H) to each Clade-II nosZ (K00376) gene.

Extends the I/II call in 05_classify_nosz.py: same DNA nhmmer machinery and the
same 8 reference alignments, but here we keep the FULL per-sub-clade bit-score
vector for every gene and report, for each Clade-II gene, its best sub-clade, the
runner-up, and the best-minus-runner-up margin (sub-clade resolution confidence).

Sub-clade reference taxa (from METHODS_nosz_clades.md):
  II.A Gemmatimonadetes        II.B Prevotella/Gemmatimonas   II.C Bacteroidetes
  II.D Dechloromonas/Dechlorosoma  II.E Desulfomonile/Desulfitobacterium
  II.G Desulfitobacterium+env   II.H Desulfosporosinus/Desulfitobacterium (Firmicutes)

Writes data/nosz_cladeII_subclades.tsv and a full score matrix to stdout.
"""
import os, csv, glob, math
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock, MSAFile
from pyhmmer.plan7 import Builder, Background

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(WORK, 'data')
LIT  = glob.glob(os.path.join(WORK, 'literature_claude*'))[0]
FA   = os.path.join(DATA, 'combined_ref.fa')

CLADE_FILES = {
    'CladeI': 'NosZREF-CladeI.fasta.text',
    'II.A':   'NosZREF-1577A.fasta.text', 'II.B': 'NosZREF-1577B.fasta.text',
    'II.C':   'NosZREF-1577C.fasta.text', 'II.D': 'NosZREF-1577D.fasta.text',
    'II.E':   'NosZREF_1577E.fasta.text', 'II.G': 'NosZREF_1577G.fasta.text',
    'II.H':   'NosZREF-1577H.fasta.text',
}
SUBCLADES = ['II.A', 'II.B', 'II.C', 'II.D', 'II.E', 'II.G', 'II.H']

# existing I/II call (so we only assign sub-clades to Clade-II genes)
clade_call, genus = {}, {}
with open(os.path.join(DATA, 'nosz_clades.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        clade_call[r['gene_id']] = r['clade']; genus[r['gene_id']] = r['genus']

# ── nosZ loci + CDS ─────────────────────────────────────────────────────────
genes = []
with open(os.path.join(DATA, 'target_genes.bed')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['ko'] == 'K00376':
            genes.append((r['prefixed_contig'], int(r['start']), int(r['end']), r['gene_id'], r['mag']))
need = {g[0] for g in genes}
seqs, cur, keep, buf = {}, None, False, []
with open(FA) as f:
    for line in f:
        if line.startswith('>'):
            if keep: seqs[cur] = ''.join(buf)
            cur = line[1:].split()[0]; keep = cur in need; buf = []
        elif keep:
            buf.append(line.strip())
    if keep: seqs[cur] = ''.join(buf)

dna = Alphabet.dna()
queries = []
for contig, s, e, gid, mag in genes:
    sub = seqs[contig][s:e].upper().replace('U', 'T')
    queries.append(TextSequence(name=gid.encode(), sequence=sub).digitize(dna))
targets = DigitalSequenceBlock(dna, queries)

# ── build HMMs + nhmmer ─────────────────────────────────────────────────────
bg, builder, hmms = Background(dna), Builder(dna), []
for cname, fn in CLADE_FILES.items():
    with MSAFile(os.path.join(LIT, fn), format='afa', digital=True, alphabet=dna) as mf:
        msa = mf.read()
    msa.name = cname.encode()
    hmm, _, _ = builder.build_msa(msa, bg); hmms.append(hmm)

def _s(x): return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
score = {g[3]: {} for g in genes}
for cname, tophits in zip(CLADE_FILES, pyhmmer.hmmer.nhmmer(hmms, targets)):
    for hit in tophits:
        gid = _s(hit.name)
        if gid in score:
            score[gid][cname] = max(score[gid].get(cname, -math.inf), hit.score)

# ── assign sub-clade to Clade-II genes ──────────────────────────────────────
rows = []
for contig, s, e, gid, mag in genes:
    if clade_call.get(gid) != 'II':
        continue
    sub_scores = {k: score[gid].get(k, 0.0) for k in SUBCLADES}
    ranked = sorted(sub_scores.items(), key=lambda kv: kv[1], reverse=True)
    best, best_s = ranked[0]
    second, second_s = ranked[1]
    margin = best_s - second_s
    conf = 'high' if margin >= 40 else 'medium' if margin >= 10 else 'low'
    rows.append({'mag': mag, 'gene_id': gid, 'genus': genus.get(gid, ''),
                 'subclade': best, 'subclade_score': f'{best_s:.1f}',
                 'runner_up': second, 'runner_up_score': f'{second_s:.1f}',
                 'subclade_margin': f'{margin:.1f}', 'confidence': conf,
                 **{k: f'{sub_scores[k]:.1f}' for k in SUBCLADES}})

cols = ['mag', 'gene_id', 'genus', 'subclade', 'subclade_score', 'runner_up',
        'runner_up_score', 'subclade_margin', 'confidence'] + SUBCLADES
with open(os.path.join(DATA, 'nosz_cladeII_subclades.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
    w.writeheader(); w.writerows(rows)

print(f"{len(rows)} Clade-II nosZ genes — sub-clade assignment\n")
print(f"  {'MAG':16s} {'genus':16s} {'sub':5s} {'score':>7s} {'2nd':>5s} {'margin':>7s}  conf")
print("  " + "-" * 70)
for r in sorted(rows, key=lambda r: (r['subclade'], r['mag'])):
    print(f"  {r['mag']:16s} {r['genus']:16s} {r['subclade']:5s} {r['subclade_score']:>7s} "
          f"{r['runner_up']:>5s} {r['subclade_margin']:>7s}  {r['confidence']}")
from collections import Counter
c = Counter(r['subclade'] for r in rows)
print("\n  sub-clade tally:", ", ".join(f"{k}:{v}" for k, v in sorted(c.items())))
print("\n=== full per-sub-clade bit-score matrix ===")
print("  " + "MAG".ljust(16) + "".join(f"{k:>8s}" for k in SUBCLADES))
for r in sorted(rows, key=lambda r: r['mag']):
    print("  " + r['mag'].ljust(16) + "".join(f"{r[k]:>8s}" for k in SUBCLADES))
print(f"\nwrote {DATA}/nosz_cladeII_subclades.tsv")
