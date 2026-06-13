#!/usr/bin/env python3
"""
Classify each MAG's nosZ (K00376) gene as Clade I vs Clade II using the
reference alignments in literature_claude*/ (Clade I = NosZREF-CladeI;
Clade II = the atypical/non-denitrifier sub-clades NosZREF-1577 A–H).

Method (DNA-level, no frame/strand assumptions):
  1. extract the 13 nosZ CDS from data/combined_ref.fa via target_genes.bed,
  2. build one nucleotide HMM per clade reference alignment (pyhmmer),
  3. nhmmer-score each gene against every clade HMM (both strands),
  4. call Clade I if the Clade-I HMM beats the best Clade-II sub-clade HMM.

Writes data/nosz_clades.tsv and prints a summary table.
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
    'II.A':   'NosZREF-1577A.fasta.text',
    'II.B':   'NosZREF-1577B.fasta.text',
    'II.C':   'NosZREF-1577C.fasta.text',
    'II.D':   'NosZREF-1577D.fasta.text',
    'II.E':   'NosZREF_1577E.fasta.text',
    'II.G':   'NosZREF_1577G.fasta.text',
    'II.H':   'NosZREF-1577H.fasta.text',
}

# ── 1. nosZ gene loci ───────────────────────────────────────────────────────
genes = []   # (contig, start, end, gene_id, mag)
with open(os.path.join(DATA, 'target_genes.bed')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['ko'] == 'K00376':
            genes.append((r['prefixed_contig'], int(r['start']), int(r['end']),
                          r['gene_id'], r['mag']))
need = {g[0] for g in genes}
print(f"{len(genes)} nosZ (K00376) genes across {len({g[4] for g in genes})} MAGs")

# genus labels
genus = {}
with open(os.path.join(DATA, 'taxonomy_labels.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        genus[r['MAG']] = r.get('Genus', '') or r.get('mag_iterativeID', '')

# ── 2. pull needed contigs from combined_ref.fa (first-token match) ─────────
seqs, cur, keep, buf = {}, None, False, []
with open(FA) as f:
    for line in f:
        if line.startswith('>'):
            if keep:
                seqs[cur] = ''.join(buf)
            cur = line[1:].split()[0]
            keep = cur in need
            buf = []
        elif keep:
            buf.append(line.strip())
    if keep:
        seqs[cur] = ''.join(buf)
missing = need - set(seqs)
if missing:
    raise SystemExit(f"contigs not found in combined_ref.fa: {missing}")

dna = Alphabet.dna()
queries = []
for contig, s, e, gid, mag in genes:
    sub = seqs[contig][s:e].upper().replace('U', 'T')
    queries.append(TextSequence(name=gid.encode(), sequence=sub).digitize(dna))
    print(f"  {mag:16s} {gid:32s} len={e-s:5d}")
targets = DigitalSequenceBlock(dna, queries)

# ── 3. build one nucleotide HMM per clade alignment ─────────────────────────
bg = Background(dna)
builder = Builder(dna)
hmms = []
for cname, fn in CLADE_FILES.items():
    with MSAFile(os.path.join(LIT, fn), format='afa', digital=True, alphabet=dna) as mf:
        msa = mf.read()
    msa.name = cname.encode()
    hmm, _, _ = builder.build_msa(msa, bg)
    hmms.append(hmm)
    print(f"  built HMM {cname:7s} (M={hmm.M}, n={len(msa.sequences)})")

# ── 4. nhmmer each clade HMM vs the nosZ genes ──────────────────────────────
def _s(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
score = {g[3]: {} for g in genes}     # gene_id -> clade -> bitscore
clade_names = list(CLADE_FILES)       # nhmmer returns TopHits in query order
for cname, tophits in zip(clade_names, pyhmmer.hmmer.nhmmer(hmms, targets)):
    for hit in tophits:
        gid = _s(hit.name)
        if gid in score:
            score[gid][cname] = max(score[gid].get(cname, -math.inf), hit.score)

# ── 5. classify ─────────────────────────────────────────────────────────────
rows = []
for contig, s, e, gid, mag in genes:
    sc = score[gid]
    cI = sc.get('CladeI', -math.inf)
    II = {k: v for k, v in sc.items() if k != 'CladeI'}
    cII = max(II.values()) if II else -math.inf
    sub = max(II, key=II.get) if II else '-'
    call = 'I' if cI >= cII else 'II'
    margin = (cI - cII) if (cI > -math.inf and cII > -math.inf) else (
        cI if cI > -math.inf else -cII)
    rows.append({'mag': mag, 'gene_id': gid, 'genus': genus.get(mag, ''),
                 'len_bp': e - s, 'clade': call,
                 'confidence': 'high' if abs(margin) >= 40 else 'low',
                 'cladeI_score': f"{cI:.1f}" if cI > -math.inf else 'NA',
                 'cladeII_score': f"{cII:.1f}" if cII > -math.inf else 'NA',
                 'cladeII_best_subclade': sub,
                 'score_margin': f"{margin:.1f}"})

with open(os.path.join(DATA, 'nosz_clades.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)

print("\n  MAG               genus            clade   I_score  II_score(sub)   margin")
print("  " + "-" * 78)
for r in sorted(rows, key=lambda r: (r['clade'], r['mag'])):
    print(f"  {r['mag']:17s} {r['genus']:15s} {('Clade '+r['clade']):8s} "
          f"{r['cladeI_score']:>7s}  {r['cladeII_score']:>7s} ({r['cladeII_best_subclade']:>4s})  "
          f"{r['score_margin']:>7s}")
nI = sum(1 for r in rows if r['clade'] == 'I')
print(f"\n  => Clade I: {nI}   Clade II: {len(rows)-nI}")
print(f"  wrote {DATA}/nosz_clades.tsv")

# ── 6. clade-split per-sample RPKM matrix ───────────────────────────────────
import statistics
from collections import defaultdict
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']
KO_GROUPS = {'nirB/D': {'K00362', 'K00363'}, 'nirS': {'K15864'},
             'nirK': {'K00368'}, 'norC/nosZ': {'K02305', 'K00376'}}
def grp_of(ko):
    return next((g for g, ks in KO_GROUPS.items() if ko in ks), None)

raw = list(csv.DictReader(open(os.path.join(DATA, 'read_counts_raw.tsv')), delimiter='\t'))
persamp = {(r['ko_group'], r['mag']): {s: float(r[s]) for s in SAMPLES}
           for r in csv.DictReader(open(os.path.join(DATA, 'gene_rpkm_per_sample.tsv')),
                                   delimiter='\t')}

# recover per-sample total mapped reads from single-gene KO groups
by_grp_mag = defaultdict(list)
for r in raw:
    g = grp_of(r['ko'])
    if g:
        by_grp_mag[(g, r['mag'])].append(r)
tot_obs = {s: [] for s in SAMPLES}
for (g, mag), gl in by_grp_mag.items():
    if len(gl) != 1 or (g, mag) not in persamp:
        continue
    L = int(gl[0]['gene_len_bp']); rp = persamp[(g, mag)]
    for s in SAMPLES:
        c = int(gl[0][s])
        if c > 0 and rp[s] > 0:
            tot_obs[s].append(c / (L / 1000.0) / rp[s] * 1e6)
TOT = {s: statistics.median(v) for s, v in tot_obs.items() if v}
print("\nrecovered per-sample total mapped reads:")
for s in SAMPLES:
    v = tot_obs[s]
    print(f"  {s}: {TOT[s]:>12,.0f}   ({len(v)} single-gene obs, "
          f"spread {min(v):,.0f}–{max(v):,.0f})")

def rpkm(c, L, s):
    return (c / (L / 1000.0)) / (TOT[s] / 1e6) if TOT.get(s) and L else 0.0

clade_of = {r['gene_id']: r['clade'] for r in rows}
OUT_ROWS = ['nirB/D', 'nirS', 'nirK', 'norC', 'nosZ_cladeI', 'nosZ_cladeII']
def in_row(rowname, r):
    ko = r['ko']
    return ((rowname == 'nirB/D'       and ko in ('K00362', 'K00363')) or
            (rowname == 'nirS'         and ko == 'K15864') or
            (rowname == 'nirK'         and ko == 'K00368') or
            (rowname == 'norC'         and ko == 'K02305') or
            (rowname == 'nosZ_cladeI'  and ko == 'K00376' and clade_of.get(r['gene_id']) == 'I') or
            (rowname == 'nosZ_cladeII' and ko == 'K00376' and clade_of.get(r['gene_id']) == 'II'))

mag_list = [l.strip() for l in open(os.path.join(DATA, 'selected_mag_list.txt')) if l.strip()]
out = [['ko_group', 'mag'] + SAMPLES]
for rowname in OUT_ROWS:
    for mag in mag_list:
        gl = [r for r in raw if r['mag'] == mag and in_row(rowname, r)]
        vals = [max((rpkm(int(r[s]), int(r['gene_len_bp']), s) for r in gl), default=0.0)
                for s in SAMPLES]
        out.append([rowname, mag] + [f"{v:.4f}" for v in vals])
with open(os.path.join(DATA, 'gene_rpkm_per_sample_claded.tsv'), 'w', newline='') as f:
    csv.writer(f, delimiter='\t').writerows(out)

# sanity: recomputed nirB/D/nirS/nirK must match the existing group-max file
mx = 0.0
for r in out[1:]:
    g, mag = r[0], r[1]
    if g in ('nirB/D', 'nirS', 'nirK') and (g, mag) in persamp:
        for i, s in enumerate(SAMPLES):
            mx = max(mx, abs(float(r[2 + i]) - persamp[(g, mag)][s]))
print(f"\n  recompute check vs existing per-sample RPKM: max abs diff = {mx:.4f}")
print(f"  wrote {DATA}/gene_rpkm_per_sample_claded.tsv  ({len(OUT_ROWS)} rows × {len(mag_list)} MAGs)")
