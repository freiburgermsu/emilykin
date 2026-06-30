#!/usr/bin/env python3
"""
Functional (KOfam) annotation of the 7 species-novel MAGs — the genomic basis
for applying the Wang/Gao/Wells supplementary characterization (emi15486-sup-
0001-supinfo.docx) to this project's NOVEL bacteria.

The reference supinfo characterises each focal organism's *functional ecology*:
carbon-source activation (acetate/propionate), storage polymers (glycogen, PHA),
the GAO gluconeogenesis route (gapA/pckA/ppdk), polyphosphate/phosphate (PAO),
type IV pili, vitamin-B12 uptake, cofactor (biotin) biosynthesis, and secreted
polysaccharide-degrading enzymes.  This script reproduces that lens *genomically*
for the 7 novel MAGs by running KOfam (pyhmmer, adaptive per-KO thresholds) over
their Bakta proteomes for a curated marker-KO panel spanning those categories.

Inputs
  ../kofamscan/all_mags.faa          Bakta proteomes, headers '>{MAG}|{locus}'
  ../kofam_new/ko_list               full KOfam thresholds (knum, threshold, score_type)
  ../kofam_new/profiles.tar.gz       full KOfam profile HMMs (members 'profiles/Kxxxxx.hmm')
Outputs (./data)
  novel_functional_ko_copynumber.tsv  MAG x KO copy-number matrix
  novel_functional_ko_long.tsv        MAG, category, KO, gene, copies (tidy)
  novel_functional_ko_hits.tsv        per-protein hits (traceability)

NOTE: metagenomic gene *potential*, not expression (no metatranscriptome). CAZy/
MEROPS/eggNOG deep annotation was not run for this dataset, so only KO-resolvable
markers are reported; secreted-protease and EPS inventories are out of reach here.
"""
import csv, os, tarfile
import pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet
from pyhmmer.plan7 import HMMFile
from pyhmmer.hmmer import hmmsearch

WORK = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WORK)
DATA = os.path.join(WORK, 'data')
SCRATCH = '/tmp/claude-1000/-home-freiburger-Documents-EmilyKin/86193b91-454f-4c29-aed9-ea9914dffd8a/scratchpad/novel_func'
os.makedirs(SCRATCH, exist_ok=True)

FAA      = os.path.join(ROOT, 'kofamscan', 'all_mags.faa')
KO_LIST  = os.path.join(ROOT, 'kofam_new', 'ko_list')
PROF_TAR = os.path.join(ROOT, 'kofam_new', 'profiles.tar.gz')

NOVEL = ['coasm_bin.55', 'CAN_2_bin.6', 'CAN_1_bin.210', 'CAN_1_bin.77',
         'CAN_5_bin.112', 'coasm_bin.260', 'coasm_bin.481']

# ── Curated marker-KO panel (category, KO, gene) — all validated vs ko_list ──
PANEL = [
    # Carbon-source activation
    ('Acetate uptake',        'K00925', 'ackA'),
    ('Acetate uptake',        'K00625', 'pta'),
    ('Acetate uptake',        'K01895', 'acs'),
    ('Propionate uptake',     'K01908', 'prpE'),
    ('Propionate uptake',     'K01659', 'prpC/mcsA'),
    ('Propionate uptake',     'K01720', 'prpD'),
    ('Propionate uptake',     'K01847', 'mcm'),
    ('Propionate uptake',     'K01848', 'mcm-N'),
    # Storage polymers — glycogen
    ('Glycogen synthesis',    'K00975', 'glgC'),
    ('Glycogen synthesis',    'K00703', 'glgA'),
    ('Glycogen synthesis',    'K00700', 'glgB'),
    ('Glycogen degradation',  'K00688', 'glgP'),
    ('Glycogen degradation',  'K01214', 'glgX'),
    ('Glycogen degradation',  'K00705', 'malQ'),
    # Storage polymers — PHA
    ('PHA cycling',           'K00626', 'phaA'),
    ('PHA cycling',           'K00023', 'phaB'),
    ('PHA cycling',           'K03821', 'phaC'),
    ('PHA cycling',           'K05973', 'phaZ'),
    # GAO gluconeogenesis (supinfo Table S3 bolded genes)
    ('Gluconeogenesis',       'K00134', 'gapA'),
    ('Gluconeogenesis',       'K01596', 'pckA(GTP)'),
    ('Gluconeogenesis',       'K01610', 'pckA(ATP)'),
    ('Gluconeogenesis',       'K01006', 'ppdk'),
    ('Gluconeogenesis',       'K03841', 'fbp'),
    ('Gluconeogenesis',       'K01007', 'pps'),
    # Polyphosphate (PAO)
    ('Polyphosphate',         'K00937', 'ppk1'),
    ('Polyphosphate',         'K22468', 'ppk2'),
    ('Polyphosphate',         'K01514', 'ppx'),
    # Phosphate transport / regulation
    ('Phosphate transport',   'K02040', 'pstS'),
    ('Phosphate transport',   'K02037', 'pstC'),
    ('Phosphate transport',   'K02038', 'pstA'),
    ('Phosphate transport',   'K02036', 'pstB'),
    ('Phosphate transport',   'K07657', 'phoB'),
    ('Phosphate transport',   'K07636', 'phoR'),
    # Type IV pili
    ('Type IV pili',          'K02650', 'pilA'),
    ('Type IV pili',          'K02652', 'pilB'),
    ('Type IV pili',          'K02653', 'pilC'),
    ('Type IV pili',          'K02654', 'pilD'),
    ('Type IV pili',          'K02662', 'pilM'),
    ('Type IV pili',          'K02663', 'pilN'),
    ('Type IV pili',          'K02664', 'pilO'),
    ('Type IV pili',          'K02665', 'pilP'),
    ('Type IV pili',          'K02666', 'pilQ'),
    ('Type IV pili',          'K02669', 'pilT'),
    ('Type IV pili',          'K02670', 'pilU'),
    # Vitamin B12 transport
    ('B12 transport',         'K16092', 'btuB'),
    ('B12 transport',         'K06073', 'btuC'),
    ('B12 transport',         'K06074', 'btuD'),
    ('B12 transport',         'K06858', 'btuF'),
    ('B12 transport',         'K03832', 'tonB'),
    ('B12 transport',         'K03561', 'exbB'),
    ('B12 transport',         'K03559', 'exbD'),
    # Vitamin B12 (de novo) biosynthesis markers
    ('B12 biosynthesis',      'K09882', 'cobS-chel'),
    ('B12 biosynthesis',      'K02232', 'cobQ'),
    ('B12 biosynthesis',      'K02233', 'cobS'),
    ('B12 biosynthesis',      'K00768', 'cobT'),
    ('B12 biosynthesis',      'K02231', 'cobU'),
    # Biotin / cofactor biosynthesis
    ('Biotin biosynthesis',   'K01012', 'bioB'),
    ('Biotin biosynthesis',   'K00833', 'bioA'),
    ('Biotin biosynthesis',   'K01935', 'bioD'),
    ('Biotin biosynthesis',   'K00652', 'bioF'),
    # Glycoside hydrolases (secreted polysaccharide breakdown — scavenger theme)
    ('Glycoside hydrolases',  'K01176', 'amyA'),
    ('Glycoside hydrolases',  'K01179', 'celA'),
    ('Glycoside hydrolases',  'K01188', 'bglX'),
    ('Glycoside hydrolases',  'K01198', 'xynB'),
    ('Glycoside hydrolases',  'K01185', 'lyz'),
    ('Glycoside hydrolases',  'K01183', 'chiA'),
]
PANEL_KOS = [ko for _, ko, _ in PANEL]
GENE = {ko: g for _, ko, g in PANEL}

# ── 1. thresholds from ko_list ──────────────────────────────────────────────
thr, stype = {}, {}
with open(KO_LIST) as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        k = r['knum']
        if k in GENE:
            t = r['threshold'].strip()
            if t and t != '-':
                thr[k] = float(t)
                stype[k] = r['score_type'].strip()
missing = [k for k in PANEL_KOS if k not in thr]
if missing:
    print('WARNING: no usable threshold for', missing)

# ── 2. extract the needed KO HMMs from profiles.tar.gz (single pass) ─────────
need = {f'profiles/{k}.hmm' for k in thr}
got = set()
print(f'extracting {len(need)} HMMs from {os.path.basename(PROF_TAR)} ...')
with tarfile.open(PROF_TAR, 'r:gz') as tf:
    for m in tf:
        if m.name in need:
            data = tf.extractfile(m).read()
            with open(os.path.join(SCRATCH, os.path.basename(m.name)), 'wb') as out:
                out.write(data)
            got.add(m.name)
            if got == need:
                break
print(f'  extracted {len(got)}/{len(need)} HMMs')

# ── 3. subset proteome to the 7 novel MAGs ──────────────────────────────────
sub_faa = os.path.join(SCRATCH, 'novel.faa')
nseq = 0
target = set(NOVEL)
with open(FAA) as fin, open(sub_faa, 'w') as fout:
    keep = False
    for line in fin:
        if line.startswith('>'):
            keep = line[1:].split('|', 1)[0] in target
            if keep:
                nseq += 1
        if keep:
            fout.write(line)
print(f'subset proteome: {nseq} proteins from {len(target)} MAGs')

# ── 4. KOfam search with adaptive thresholds ─────────────────────────────────
alpha = Alphabet.amino()
with SequenceFile(sub_faa, digital=True, alphabet=alpha) as sf:
    seqs = sf.read_block()

hmms = []
for k in thr:
    with HMMFile(os.path.join(SCRATCH, f'{k}.hmm')) as hf:
        hm = hf.read()
        hmms.append(hm)

# map HMM NAME -> KO (KOfam HMM NAME == KO id; verify & fall back to accession)
def as_str(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else (x or '')

copy = {m: {k: 0 for k in thr} for m in NOVEL}
hit_rows = []
for hits in hmmsearch(hmms, seqs, cpus=os.cpu_count() or 4):
    ko = as_str(hits.query.name)
    if ko not in thr:
        # fallback: some profiles carry the KO as the accession
        acc = as_str(hits.query.accession)
        ko = acc if acc in thr else ko
    if ko not in thr:
        continue
    t, st = thr[ko], stype[ko]
    for h in hits:
        score = h.score if st == 'full' else (h.best_domain.score if h.best_domain else h.score)
        if score >= t:
            name = as_str(h.name)
            mag = name.split('|', 1)[0]
            if mag in copy:
                copy[mag][ko] += 1
                hit_rows.append((mag, ko, GENE[ko], name, f'{score:.1f}', st, f'{t:.2f}'))

# ── 5. write outputs ─────────────────────────────────────────────────────────
os.makedirs(DATA, exist_ok=True)
kos = [k for k in PANEL_KOS if k in thr]

with open(os.path.join(DATA, 'novel_functional_ko_copynumber.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['MAG'] + kos)
    for m in NOVEL:
        w.writerow([m] + [copy[m][k] for k in kos])

with open(os.path.join(DATA, 'novel_functional_ko_long.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['MAG', 'category', 'KO', 'gene', 'copies'])
    for cat, ko, g in PANEL:
        if ko in thr:
            for m in NOVEL:
                w.writerow([m, cat, ko, g, copy[m][ko]])

with open(os.path.join(DATA, 'novel_functional_ko_hits.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['MAG', 'KO', 'gene', 'protein_id', 'score', 'score_type', 'threshold'])
    for row in sorted(hit_rows):
        w.writerow(row)

# ── 6. console summary by category ───────────────────────────────────────────
print('\n=== copy-number summary (rows=MAG) ===')
cats = []
for cat, _, _ in PANEL:
    if cat not in cats:
        cats.append(cat)
for m in NOVEL:
    present = {cat: 0 for cat in cats}
    for cat, ko, g in PANEL:
        if ko in thr and copy[m][ko] > 0:
            present[cat] += 1
    line = '  '.join(f'{cat}:{present[cat]}' for cat in cats if present[cat])
    print(f'{m:16s} {line}')
print(f'\nwrote 3 TSVs to {DATA}')
