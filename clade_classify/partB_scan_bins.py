#!/usr/bin/env python3
"""
Part B: Scan all 276 bins for clade III L-NosZ (He et al. 2025 N2OR).

Uses the He 269NosZ HMM (built in Part A) to search all bin nucleotide contigs.
Translated hits are screened with CuA / CuZ motif:
  CuA = C..FC...H.EM  (present in all N2ORs)
  CuZ = G.HH (L-NosZ / clade III) vs D.HH (C-NosZ / clade I-II)
Confirmed L-NosZ candidates also checked by diamond vs Chee+Orellana (<35% id).

Output: clade_classify/out/partB_motif.tsv, partB_lnosz_hits.tsv
"""
import csv, re, subprocess
from pathlib import Path

HERE    = Path(__file__).parent
OUT     = HERE / 'out'; OUT.mkdir(exist_ok=True)
HMM     = OUT / 'He_269NosZ.hmm'
DREP    = Path('/scratch1/afreiburger/emilykin/processed/mag/drep/dereplicated_genomes')
CHEE    = HERE / 'Chee_plus_Orellana_NosZ.prot.raw.faa'

ENV2    = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_/bin'
HMMSEARCH = f'{ENV2}/hmmsearch'
DIAMOND   = f'{ENV2}/diamond'

PYBIN   = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin/python3'

CUA  = re.compile(r'C.{2}FC.{3}H.EM')
CUZC = re.compile(r'D.HH')   # C-NosZ
CUZL = re.compile(r'G.HH')   # L-NosZ (clade III)

def read_fasta(path):
    seqs = {}; cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                cur = line[1:].split()[0]; seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}

def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for name, seq in seqs.items():
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def translate_dna(dna):
    codon = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    aas = []
    for strand, seq in [('+', dna), ('-', dna[::-1].translate(comp))]:
        for frame in range(3):
            s = seq[frame:]; parts = s[:len(s)-(len(s)%3)]
            aa = ''.join(codon.get(parts[i:i+3].upper(),'X') for i in range(0,len(parts),3))
            for j, frag in enumerate(aa.split('*')):
                if len(frag) >= 100:
                    aas.append((f'fr{frame}{strand}orf{j}', frag))
    return aas

if not HMM.exists():
    print(f'ERROR: {HMM} not found. Run classify_nosz.py (Part A) first.')
    raise SystemExit(1)

# ── B1. Build concatenated bin FASTA (if not already done) ───────────────────
all_bins = OUT / 'all_bins.fna'
bin_fas  = sorted(DREP.glob('*.fa'))
print(f'Found {len(bin_fas)} bins')

if not all_bins.exists():
    print('Building all_bins.fna...')
    with open(all_bins, 'w') as out_f:
        for fa in bin_fas:
            bn = fa.stem
            with open(fa) as inp:
                for line in inp:
                    if line.startswith('>'):
                        out_f.write(f'>{bn}::{line[1:]}')
                    else:
                        out_f.write(line)
    print(f'  Written {all_bins}')

# ── B2. hmmsearch all bins (6-frame translated) ───────────────────────────────
tblout = OUT / 'partB_hmmsearch.tblout'
if not tblout.exists():
    print('Running hmmsearch (6-frame DNA translation)...')
    subprocess.run([HMMSEARCH, '--cpu', '16', '--tblout', str(tblout),
                    '--noali', '--dnax', '-E', '1e-10', str(HMM), str(all_bins)],
                   check=True, capture_output=True)
    print(f'  Written {tblout}')
else:
    print(f'Using existing {tblout}')

# Parse hit contigs
hit_contigs = set()
hit_scores  = {}
with open(tblout) as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.split()
        if len(p) < 10: continue
        # --dnax target name: contig/frame
        name = p[0]
        score = float(p[5])
        # extract contig name (everything before the last /N frame indicator)
        contig = name.rsplit('/', 1)[0] if '/' in name else name
        hit_contigs.add(contig)
        hit_scores[name] = score

print(f'  Raw hits: {len(hit_scores)}, unique contigs: {len(hit_contigs)}')

if not hit_contigs:
    print('\nNo hits found — clade III L-NosZ not detected in bins.')
    with open(OUT / 'partB_motif.tsv', 'w') as f:
        f.write('orf\tbin\tCuA\tCuZ\tcall\thmm_score\n')
    with open(OUT / 'partB_lnosz_hits.tsv', 'w') as f:
        f.write('orf\tbin\tCuA\tCuZ\tcall\thmm_score\tpct_id_to_CNosZ\tconfirmed_LNosZ\n')
    print('Written empty output files.')
    raise SystemExit(0)

# ── B3. Extract hit contigs and get ORFs ─────────────────────────────────────
print('\nExtracting hit contigs...')
hit_fna = OUT / 'partB_hit_contigs.fna'
hit_seqs = {}
with open(all_bins) as f:
    cur = None
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            cur = line[1:].split()[0]
        elif cur and cur in hit_contigs:
            hit_seqs.setdefault(cur, []).append(line)
hit_seqs = {k: ''.join(v) for k, v in hit_seqs.items()}
write_fasta(hit_fna, hit_seqs)
print(f'  {len(hit_seqs)} hit contigs')

# ── B4. 6-frame translate + confirm with protein hmmsearch ───────────────────
print('\nTranslating hit contigs (6-frame)...')
all_orfs = {}
for cname, dna in hit_seqs.items():
    for orf_suffix, aa in translate_dna(dna):
        all_orfs[f'{cname}__{orf_suffix}'] = aa

print(f'  {len(all_orfs)} candidate ORFs')
orf_faa = OUT / 'partB_orfs.faa'
write_fasta(orf_faa, all_orfs)

print('\nConfirming NosZ ORFs with protein hmmsearch...')
conf_tblout = OUT / 'partB_confirmed.tblout'
subprocess.run([HMMSEARCH, '--cpu', '8', '--tblout', str(conf_tblout),
                '--noali', '-E', '1e-5', str(HMM), str(orf_faa)],
               check=True, capture_output=True)

confirmed = {}
with open(conf_tblout) as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.split()
        if len(p) < 10: continue
        confirmed[p[0]] = float(p[5])

print(f'  Confirmed NosZ ORFs: {len(confirmed)}')

# ── B5. CuA / CuZ motif check ─────────────────────────────────────────────────
print('\nMotif diagnosis...')
motif_rows = []
for orf_id, aa in all_orfs.items():
    if orf_id not in confirmed: continue
    bin_name = orf_id.split('::')[0] if '::' in orf_id else orf_id.split('__')[0]
    m = CUA.search(aa)
    head = aa[:m.start()] if m else aa  # CuZ is N-terminal of CuA
    c_hit = bool(CUZC.search(head))
    l_hit = bool(CUZL.search(head))
    if l_hit and not c_hit:
        call = 'L-NosZ_cladeIII'
    elif c_hit:
        call = 'C-NosZ_cladeI_II'
    else:
        # fallback: search tail too
        tail = aa[m.end():] if m else aa
        if bool(CUZC.search(tail)): call = 'C-NosZ_cladeI_II'
        elif bool(CUZL.search(tail)): call = 'L-NosZ_cladeIII'
        elif m: call = 'CuA_present_no_CuZ'
        else: call = 'no_CuA'
    cuz_str = 'GXHH' if call.startswith('L-NosZ') else ('DXHH' if 'cladeI_II' in call else '?')
    motif_rows.append({
        'orf': orf_id, 'bin': bin_name,
        'CuA': 'CuA' if m else 'noCuA', 'CuZ': cuz_str,
        'call': call, 'hmm_score': confirmed[orf_id]
    })

with open(OUT / 'partB_motif.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['orf','bin','CuA','CuZ','call','hmm_score'], delimiter='\t')
    w.writeheader(); w.writerows(motif_rows)
print(f'  Written {OUT}/partB_motif.tsv ({len(motif_rows)} NosZ ORFs)')

lnosz_cands = [r for r in motif_rows if r['call'] == 'L-NosZ_cladeIII']
print(f'\nL-NosZ motif candidates (GXHH, not DXHH): {len(lnosz_cands)}')
for r in lnosz_cands:
    print(f'  {r["orf"]}  score={r["hmm_score"]:.1f}')

# ── B6. Confirm L-NosZ: diamond vs Chee+Orellana (<35% = confirmed L-NosZ) ──
lnosz_confirmed = []
if lnosz_cands:
    print('\nDiamond check: L-NosZ candidates vs Chee+Orellana C-NosZ...')
    cands_faa = OUT / 'partB_lnosz_cands.faa'
    write_fasta(cands_faa, {r['orf']: all_orfs[r['orf']] for r in lnosz_cands})
    chee_db = OUT / 'chee_diamond.dmnd'
    if not chee_db.exists():
        subprocess.run([DIAMOND, 'makedb', '--in', str(CHEE), '--db', str(chee_db), '--quiet'],
                       check=True, capture_output=True)
    dia_out = OUT / 'partB_lnosz_vs_cnosz.tsv'
    subprocess.run([DIAMOND, 'blastp', '-q', str(cands_faa), '--db', str(chee_db),
                    '--out', str(dia_out), '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
                    '--max-target-seqs', '1', '--quiet', '--evalue', '0.001'],
                   check=True, capture_output=True)
    pcts = {}
    with open(dia_out) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 3: pcts[p[0]] = float(p[2])
    for r in lnosz_cands:
        pct = pcts.get(r['orf'], None)
        confirmed_flag = 'YES' if (pct is None or pct < 35.0) else f'NO (pct_id={pct:.1f}%)'
        lnosz_confirmed.append({**r,
            'pct_id_to_CNosZ': f'{pct:.1f}%' if pct else 'no_diamond_hit',
            'confirmed_LNosZ': confirmed_flag})

with open(OUT / 'partB_lnosz_hits.tsv', 'w', newline='') as f:
    if lnosz_confirmed:
        w = csv.DictWriter(f, fieldnames=list(lnosz_confirmed[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(lnosz_confirmed)
    else:
        f.write('orf\tbin\tCuA\tCuZ\tcall\thmm_score\tpct_id_to_CNosZ\tconfirmed_LNosZ\n')

# ── Summary ───────────────────────────────────────────────────────────────────
print('\n=== PART B SUMMARY ===')
print(f'Bins scanned:                {len(bin_fas)}')
print(f'Raw hmmsearch hits:          {len(hit_scores)}')
print(f'Confirmed NosZ ORFs:         {len(confirmed)}')
lnosz_yes = [r for r in lnosz_confirmed if r.get('confirmed_LNosZ','').startswith('YES')]
print(f'L-NosZ motif (GXHH):        {len(lnosz_cands)}')
print(f'Confirmed L-NosZ (<35% id): {len(lnosz_yes)}')
if lnosz_yes:
    print('  CLADE III L-NosZ PRESENT:')
    for r in lnosz_yes:
        print(f'    {r["bin"]}: {r["orf"]}  score={r["hmm_score"]:.1f}')
else:
    print('  Clade III L-NosZ: NOT DETECTED in any bin.')
