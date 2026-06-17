#!/usr/bin/env python3
"""
Part A (fixed): reclassify 13 K00376 nosZ genes as clade I / II.
Uses diamond similarity to Chee+Orellana reference set + CuA/CuZ motif check.

Fixes from classify_nosz.py:
  1. nosz_clades.tsv taken from repo path, not working dir.
  2. CuZ motif searched N-terminal of CuA (not C-terminal).
  3. Clade I/II from nearest Chee+Orellana hit (diamond), not dendropy PDM
     (which failed on the large alignment).
"""
import csv, json, os, re, subprocess
from pathlib import Path

HERE   = Path(__file__).parent
OUT    = HERE / 'out'; OUT.mkdir(exist_ok=True)
CHEE   = HERE / 'Chee_plus_Orellana_NosZ.prot.raw.faa'
REPO   = Path('/scratch1/afreiburger/emilykin/EmilyKin')

# The working gene_ab_figure data (may differ from repo copy)
WDATA  = Path('/scratch1/afreiburger/emilykin/gene_ab_figure/data')
# nosz_clades.tsv is only in the repo
PRIOR_TSV = REPO / 'gene_ab_figure/data/nosz_clades.tsv'

PROC   = Path('/scratch1/afreiburger/emilykin/processed/mag')
ENV2   = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_/bin'
DIAMOND = f'{ENV2}/diamond'

# ── Clade I vs II organism name lookup ───────────────────────────────────────
# Clade I = Pseudomonadota "typical denitrifier" nosZ lineage
CLADE_I_GENERA = {
    'Dechloromonas', 'Rhodoferax', 'Burkholderia', 'Thiobacillus',
    'Aromatoleum', 'Azoarcus', 'Lautropia', 'Alicycliphilus', 'Acidovorax',
    'Comamonas', 'Paracoccus', 'Pseudomonas', 'Shewanella', 'Bradyrhizobium',
    'Thauera', 'Magnetospirillum', 'Hydrogenophaga', 'Rubrivivax',
    'Sulfurihydrogenibium', 'Aquifex', 'Persephonella',
}
# Clade II = atypical (Bacteroidota, Myxococcota, Accumulibacter IIA-G, etc.)
CLADE_II_GENERA = {
    'Wolinella', 'Anaeromyxobacter', 'Gemmatimonas', 'Flavobacterium',
    'Cytophaga', 'Sphingobacterium', 'Ignavibacterium', 'Chlorobaculum',
    'Accumulibacter', 'Denitratisoma', 'Azonexus', 'Azospira',
    'Magnetococcus', 'Rhodothermus', 'Persicobacter', 'Microscilla',
    'Aureispira', 'Gramella', 'Kordia',
}

def classify_hit_name(name):
    """Return 'I', 'II', or None from the Chee+Orellana hit name."""
    for genus in CLADE_I_GENERA:
        if genus.lower() in name.lower():
            return 'I'
    for genus in CLADE_II_GENERA:
        if genus.lower() in name.lower():
            return 'II'
    return None

# ── Motif patterns ────────────────────────────────────────────────────────────
CUA  = re.compile(r'C.{2}FC.{3}H.EM')  # CuA — all N2ORs
CUZC = re.compile(r'D.HH')             # CuZ — C-NosZ (clade I/II)
CUZL = re.compile(r'G.HH')             # CuZ — L-NosZ (clade III)

# ── Read pre-built query proteins ─────────────────────────────────────────────
query_faa = OUT / 'nosz13_query.faa'
if not query_faa.exists():
    print(f'ERROR: {query_faa} not found — run classify_nosz.py first to generate it.')
    raise SystemExit(1)

query_prots = {}
with open(query_faa) as f:
    cur = None
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            cur = line[1:].split()[0]; query_prots[cur] = []
        elif cur:
            query_prots[cur].append(line)
query_prots = {k: ''.join(v) for k, v in query_prots.items()}
print(f'Loaded {len(query_prots)} query nosZ proteins from {query_faa}')

# safe_id → original gene_id mapping
# safe_id = gid.replace('|','_').replace('.','_')
def safe_to_gid(safe_id, all_gids):
    for g in all_gids:
        if g.replace('|','_').replace('.','_') == safe_id:
            return g
    return safe_id

# ── Load K00376 gene metadata ─────────────────────────────────────────────────
nosz_meta = {}  # safe_id -> {gene_id, mag}
with open(WDATA / 'target_genes.bed') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['ko'] == 'K00376':
            gid   = row['gene_id']
            safe  = gid.replace('|','_').replace('.','_')
            nosz_meta[safe] = {'gene_id': gid, 'mag': row['mag']}

# ── Load prior nosz_clades.tsv ────────────────────────────────────────────────
prior = {}
if PRIOR_TSV.exists():
    with open(PRIOR_TSV) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            prior[row['gene_id']] = row
    print(f'Loaded {len(prior)} prior clade calls from {PRIOR_TSV}')
else:
    print(f'WARNING: {PRIOR_TSV} not found')

# ── A1. diamond blastp: query vs Chee+Orellana ───────────────────────────────
print('\nA1. Running diamond blastp: query nosZ vs Chee+Orellana reference...')
chee_db = OUT / 'chee_diamond.dmnd'
if not chee_db.exists():
    subprocess.run([DIAMOND, 'makedb', '--in', str(CHEE), '--db', str(chee_db), '--quiet'],
                   check=True, capture_output=True)
    print(f'  Built diamond db: {chee_db}')

diamond_out = OUT / 'nosz13_vs_chee.tsv'
subprocess.run([
    DIAMOND, 'blastp', '-q', str(query_faa), '--db', str(chee_db),
    '--out', str(diamond_out), '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue',
    '--max-target-seqs', '5', '--quiet', '--evalue', '1e-10',
], check=True, capture_output=True)

# Parse top hits per query
diamond_hits = {}  # safe_id -> [(sseqid, pident, length, evalue)]
with open(diamond_out) as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p) < 5: continue
        q = p[0]
        diamond_hits.setdefault(q, []).append((p[1], float(p[2]), int(p[3]), float(p[4])))

print(f'  Queries with hits: {len(diamond_hits)}')

# ── A2. Assign clade from top hits ───────────────────────────────────────────
print('\nA2. Assigning clade I / II from diamond nearest-neighbor...')
diamond_clade = {}  # safe_id -> ('I'/'II'/'ambiguous', top_hit, pident)
for safe, hits in diamond_hits.items():
    calls = []
    for sseqid, pident, length, evalue in hits[:5]:
        c = classify_hit_name(sseqid)
        if c:
            calls.append(c)
    if calls:
        # majority vote
        c1 = calls.count('I'); c2 = calls.count('II')
        if c1 > c2: clade = 'I'
        elif c2 > c1: clade = 'II'
        else: clade = 'ambiguous'
    else:
        clade = 'unknown'
    top = hits[0]
    diamond_clade[safe] = (clade, top[0], top[1])

# ── A3. CuA / CuZ motif check ────────────────────────────────────────────────
print('\nA3. CuA / CuZ motif check...')
motif_result = {}  # safe_id -> dict
for safe, aa in query_prots.items():
    m = CUA.search(aa)
    if m:
        # CuZ is N-terminal of CuA in NosZ — search the region BEFORE CuA
        head = aa[:m.start()]
    else:
        head = aa  # fallback: search whole sequence
    c_hit = bool(CUZC.search(head))
    l_hit = bool(CUZL.search(head))
    if l_hit and not c_hit:
        call = 'L-NosZ_cladeIII'
    elif c_hit:
        call = 'C-NosZ_cladeI_II'
    else:
        # Try searching C-terminal of CuA as fallback (some variants)
        tail = aa[m.end():] if m else aa
        c2 = bool(CUZC.search(tail)); l2 = bool(CUZL.search(tail))
        if l2 and not c2:
            call = 'L-NosZ_cladeIII'
        elif c2:
            call = 'C-NosZ_cladeI_II (C-term)'
        elif m:
            call = 'CuA_present_no_CuZ'
        else:
            call = 'no_CuA'
    motif_result[safe] = {
        'CuA': 'CuA' if m else 'noCuA',
        'CuZ': 'GXHH' if (call.startswith('L-NosZ')) else
               'DXHH' if 'cladeI_II' in call else '?',
        'call': call,
    }
    print(f'  {safe}: {motif_result[safe]["CuA"]}  CuZ={motif_result[safe]["CuZ"]}  → {call}')

# ── A4. Consensus clade + comparison with prior ───────────────────────────────
print('\nA4. Consensus and comparison with prior calls...')
rows = []
changed = []
for safe, info in nosz_meta.items():
    gid  = info['gene_id']
    mag  = info['mag']
    pr   = prior.get(gid, {})
    prior_clade = pr.get('clade', 'N/A')

    dc   = diamond_clade.get(safe, ('N/A', 'N/A', 0.0))
    dia_clade, top_hit, pident = dc
    mot  = motif_result.get(safe, {})

    # Motif is deterministic for clade III; use diamond for I/II
    if mot.get('call', '').startswith('L-NosZ'):
        final_clade = 'III'
    elif dia_clade in ('I', 'II'):
        final_clade = dia_clade
    else:
        final_clade = prior_clade  # fall back

    if prior_clade != 'N/A' and prior_clade != final_clade:
        changed.append((gid, prior_clade, final_clade))

    rows.append({
        'gene_id': gid, 'mag': mag,
        'prior_clade': prior_clade,
        'diamond_clade': dia_clade,
        'top_chee_hit': top_hit,
        'top_pident': f'{pident:.1f}',
        'final_clade': final_clade,
        'CuA': mot.get('CuA', 'N/A'),
        'CuZ': mot.get('CuZ', 'N/A'),
        'motif_call': mot.get('call', 'N/A'),
    })

# ── Write output ──────────────────────────────────────────────────────────────
fields = ['gene_id','mag','prior_clade','diamond_clade','top_chee_hit',
          'top_pident','final_clade','CuA','CuZ','motif_call']
with open(OUT / 'partA_clade_I_II.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(f'\nWritten: {OUT}/partA_clade_I_II.tsv')

# ── Summary ───────────────────────────────────────────────────────────────────
print('\n=== PART A SUMMARY ===')
for r in rows:
    chg = ' ← CHANGED' if r['prior_clade'] != r['final_clade'] and r['prior_clade'] != 'N/A' else ''
    print(f"  {r['mag']:20s}  prior={r['prior_clade']}  diamond={r['diamond_clade']}  "
          f"final={r['final_clade']}  {r['CuZ']}  {chg}")

if changed:
    print(f'\nCHANGED from prior: {len(changed)}')
    for gid, pc, fc in changed:
        print(f'  {gid}: {pc} → {fc}')
else:
    print('\nNo changes from prior clade I/II calls.')

# ── Save updated nosz_clades.tsv if changed ───────────────────────────────────
if changed:
    updated = {}
    for r in rows:
        updated[r['gene_id']] = r['final_clade']
    # Rebuild nosz_clades.tsv with updated clade column
    prior_rows = []
    with open(PRIOR_TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            if row['gene_id'] in updated:
                row['clade'] = updated[row['gene_id']]
            prior_rows.append(row)
    updated_tsv = OUT / 'nosz_clades_updated.tsv'
    with open(updated_tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader(); w.writerows(prior_rows)
    print(f'\nUpdated nosz_clades.tsv written to: {updated_tsv}')
    print('  Copy to gene_ab_figure/data/nosz_clades.tsv if changes are accepted.')
else:
    print(f'\nExisting nosz_clades.tsv is correct — no update needed.')
