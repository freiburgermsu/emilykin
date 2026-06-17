#!/usr/bin/env python3
"""
nosZ clade I / II / III classification (Part A + Part B bins).

Part A — Reclassify the 13 K00376 C-NosZ genes as clade I vs II using:
  1. Protein extraction from combined_ref.fa (translated CDS)
  2. hmmbuild from He 269NosZ reference alignment
  3. hmmalign to profile-align the 13 queries + Chee+Orellana references
  4. FastTree to build a reference+query tree
  5. Clade assignment by nearest-neighbor reference taxa

Part B — Screen all 276 MAG bins for clade III L-NosZ:
  1. hmmsearch of all bins (nucleotide) with He HMM → hit ORFs
  2. CuA / CuZ motif diagnosis on hit ORFs
  3. Confirm L-NosZ by < 35% identity to C-NosZ (diamond)

Outputs → clade_classify/out/
  partA_clade_I_II.tsv
  partB_motif.tsv
  partB_lnosz_hits.tsv (if any)
  summary.md
"""
import csv, json, os, re, subprocess, sys, textwrap
from pathlib import Path
from collections import defaultdict

# ── Paths ────────────────────────────────────────────────────────────────────
HERE     = Path(__file__).parent
OUT      = HERE / 'out'; OUT.mkdir(exist_ok=True)
GPKG_X   = HERE / 'gpkg_x'
REF_ALN  = GPKG_X / 'graftM_269NosZ_deduplicated_aligned.fasta'
SEQINFO  = GPKG_X / 'graftM_269NosZ_seqinfo.csv'
CHEE_FAA = HERE / 'Chee_plus_Orellana_NosZ.prot.raw.faa'

WORK      = Path('/scratch1/afreiburger/emilykin/gene_ab_figure')
DATA      = WORK / 'data'
PROC      = Path('/scratch1/afreiburger/emilykin/processed/mag')
DREP      = Path('/scratch1/afreiburger/emilykin/processed/mag/drep/dereplicated_genomes')

ENV2      = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/8c703f8abdc13cdc1ef374a2bcf37f14_/bin'
HMMBUILD  = f'{ENV2}/hmmbuild'
HMMALIGN  = f'{ENV2}/hmmalign'
HMMSEARCH = f'{ENV2}/hmmsearch'
FASTTREE  = f'{ENV2}/FastTreeMP'
PPLACER   = f'{ENV2}/pplacer'
DIAMOND   = f'{ENV2}/diamond'
PYBIN     = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin/python3'

SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

# ── Helper ───────────────────────────────────────────────────────────────────
def run(cmd, **kwargs):
    print(f'  $ {" ".join(str(x) for x in cmd)}')
    r = subprocess.run([str(x) for x in cmd], check=True,
                       capture_output=True, text=True, **kwargs)
    if r.stdout.strip():
        print(r.stdout[-2000:])
    return r

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

def translate_dna(dna, strand='+'):
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
    if strand == '-':
        dna = dna[::-1].translate(comp)
    aa = []
    for i in range(0, len(dna) - 2, 3):
        aa.append(codon.get(dna[i:i+3].upper(), 'X'))
    prot = ''.join(aa)
    if prot.endswith('*'):
        prot = prot[:-1]
    return prot

# ═══════════════════════════════════════════════════════════════════════════════
# PART A — Extract and classify the 13 K00376 nosZ proteins
# ═══════════════════════════════════════════════════════════════════════════════
print('=' * 70)
print('PART A — nosZ clade I / II classification')
print('=' * 70)

# A1. Read target_genes.bed for K00376 entries
print('\nA1. Loading K00376 gene coordinates from target_genes.bed...')
nosz_genes = []
with open(DATA / 'target_genes.bed') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['ko'] == 'K00376':
            nosz_genes.append(row)
print(f'  Found {len(nosz_genes)} K00376 (nosZ) genes')

# A2. Extract protein sequences from combined_ref.fa (translate CDS)
print('\nA2. Extracting and translating nosZ CDS...')
ref_seqs = read_fasta(DATA / 'combined_ref.fa')
query_prots = {}  # gene_id -> aa_seq
for g in nosz_genes:
    prefix = g['prefixed_contig']
    s0, e  = int(g['start']), int(g['end'])
    mag    = g['mag']
    gid    = g['gene_id']
    # Determine strand from bakta GFF3
    prot_id = gid.split('|', 1)[1] if '|' in gid else gid
    strand = '+'
    gff = PROC / 'bakta' / mag / f'{mag}.gff3'
    if gff.exists():
        with open(gff) as f:
            for line in f:
                if line.startswith('#'): continue
                cols = line.split('\t')
                if len(cols) < 9 or cols[2] != 'CDS': continue
                if f'locus_tag={prot_id}' in cols[8]:
                    strand = cols[6]; break
    dna = ref_seqs.get(prefix, '')[s0:e]
    if not dna:
        print(f'  WARNING: no sequence for {prefix}')
        continue
    aa = translate_dna(dna, strand)
    safe_id = gid.replace('|', '_').replace('.', '_')
    query_prots[safe_id] = aa
    print(f'  {gid}  len={len(aa)} aa  strand={strand}')

qfaa = OUT / 'nosz13_query.faa'
write_fasta(qfaa, query_prots)
print(f'  Written: {qfaa}')

# A3. Build HMM from He 269NosZ reference alignment
print('\nA3. Building HMM from He 269NosZ reference alignment...')
hmm = OUT / 'He_269NosZ.hmm'
run([HMMBUILD, '--amino', hmm, REF_ALN])

# A4. hmmalign query proteins + Chee+Orellana to HMM
print('\nA4. hmmalign: aligning query proteins to He HMM...')
# Merge query proteins with Chee+Orellana for the tree
all_query_faa = OUT / 'all_query_plus_chee.faa'
chee_seqs = read_fasta(CHEE_FAA)
# Rename Chee sequences for clean tree labels
chee_renamed = {}
for k, v in chee_seqs.items():
    # Strip trailing newlines and spaces in sequence
    clean_k = k.split()[0].replace('|', '_').replace(' ', '_')[:60]
    chee_renamed[f'REF__{clean_k}'] = v
combined_query = {**query_prots, **chee_renamed}
write_fasta(all_query_faa, combined_query)

# Also include He reference sequences (de-gapped) for the tree
he_seqs_raw = read_fasta(REF_ALN)
he_degapped = {}
for k, v in he_seqs_raw.items():
    clean_k = k.split()[0].replace('|', '_').replace(' ', '_')[:60]
    he_degapped[f'HE__{clean_k}'] = v.replace('-', '').replace('.', '')
he_faa = OUT / 'he_refs_degapped.faa'
write_fasta(he_faa, he_degapped)

# Align all (query + Chee + He references) to the HMM
all_to_align = OUT / 'all_for_alignment.faa'
all_seqs_combined = {**query_prots, **chee_renamed, **he_degapped}
write_fasta(all_to_align, all_seqs_combined)
aln_out = OUT / 'nosz_aligned.afa'
run([HMMALIGN, '--amino', '--trim', '--outformat', 'afa', '-o', aln_out, hmm, all_to_align])

# A5. Build phylogenetic tree
print('\nA5. Building phylogenetic tree with FastTreeMP...')
tree_out = OUT / 'nosz_tree.nwk'
run([FASTTREE, '-wag', '-gamma', '-log', OUT / 'fasttree.log',
     '-out', tree_out, aln_out])

print(f'  Tree: {tree_out}')

# A6. Parse tree for clade assignments
print('\nA6. Parsing tree for clade I / II / III assignment...')

# Known clade I marker taxa (Pseudomonadota C-NosZ typical):
CLADE_I_MARKERS = [
    'Pseudomonas', 'Paracoccus', 'Bradyrhizobium', 'Azoarcus', 'Thauera',
    'Shewanella', 'Aromatoleum', 'Rhodoferax', 'Thiobacillus', 'Dechloromonas',
    'Burkholderia', 'Alicycliphilus', 'Acidovorax', 'Lautropia', 'Comamonas',
]
# Known clade II marker taxa:
CLADE_II_MARKERS = [
    'Wolinella', 'Anaeromyxobacter', 'Bacteroidota', 'Gemmatimonadota',
    'Accumulibacter', 'Flavobacteriia', 'Bacteroidia', 'Sphingobacteriia',
    'Cytophagales', 'Ignavibacteria', 'Chlorobiales',
]
# Clade III (L-NosZ) marker taxa:
CLADE_III_MARKERS = [
    'Desulfobacterota', 'Nitrospinota', 'Chloroflexota', 'Moorella',
    'Desulfitobacterium', 'Desulfosporosinus', 'Desulforamulus',
    'Thermacetogenium', 'Geobacillus', 'Bacillota_LNosZ', 'Frigididesulfovibrio',
]

# Use dendropy for tree parsing if available, else simple Newick neighbor scan
try:
    import dendropy
    HAS_DENDROPY = True
except ImportError:
    HAS_DENDROPY = False

def classify_by_nearest_ref(tree_str, query_ids, clade_i_markers, clade_ii_markers, clade_iii_markers):
    """
    Simple nearest-neighbor: for each query leaf, find its nearest labeled
    reference leaf and assign clade.
    Returns {query_id: {'clade': 'I'|'II'|'III'|'?', 'nearest_ref': str, 'dist': float}}
    """
    if not HAS_DENDROPY:
        return {}
    import dendropy
    tree = dendropy.Tree.get(data=tree_str, schema='newick')
    # Build distance matrix
    pdm = tree.phylogenetic_distance_matrix()
    results = {}
    for leaf in tree.taxon_namespace:
        name = leaf.label
        if name is None: continue
        is_query = any(qid == name for qid in query_ids)
        if not is_query: continue
        best_dist, best_ref, best_clade = float('inf'), None, '?'
        for ref in tree.taxon_namespace:
            rname = ref.label
            if rname is None or rname == name: continue
            if not (rname.startswith('REF__') or rname.startswith('HE__')):
                continue
            try:
                d = pdm(leaf, ref)
            except Exception:
                continue
            if d < best_dist:
                best_dist = d; best_ref = rname
                # Classify based on clade markers in the ref name
                rn_lower = rname.lower()
                if any(m.lower() in rn_lower for m in clade_iii_markers):
                    best_clade = 'III'
                elif any(m.lower() in rn_lower for m in clade_i_markers):
                    best_clade = 'I'
                elif any(m.lower() in rn_lower for m in clade_ii_markers):
                    best_clade = 'II'
                else:
                    best_clade = 'unknown_ref'
        results[name] = {'clade': best_clade, 'nearest_ref': best_ref, 'dist': best_dist}
    return results

if tree_out.exists():
    tree_str = tree_out.read_text()
    placements = classify_by_nearest_ref(
        tree_str, list(query_prots.keys()),
        CLADE_I_MARKERS, CLADE_II_MARKERS, CLADE_III_MARKERS)
else:
    placements = {}
    print('  WARNING: tree not generated, skipping placement')

# A7. Motif check (CuA + CuZ)
print('\nA7. Running CuA / CuZ motif diagnosis...')
CUA  = re.compile(r'C.{2}FC.{3}H.EM')
CUZC = re.compile(r'D.HH')   # C-NosZ clade I/II
CUZL = re.compile(r'G.HH')   # L-NosZ clade III

motif_results = {}
for sid, aa in query_prots.items():
    m = CUA.search(aa)
    tail = aa[m.end():] if m else aa
    c_hit = bool(CUZC.search(tail))
    l_hit = bool(CUZL.search(tail))
    if l_hit and not c_hit:
        call = 'L-NosZ_cladeIII'
    elif c_hit:
        call = 'C-NosZ_cladeI_II'
    else:
        call = 'ambiguous'
    motif_results[sid] = {
        'CuA': 'CuA' if m else 'noCuA',
        'CuZ': 'GXHH' if (l_hit and not c_hit) else ('DXHH' if c_hit else '?'),
        'motif_call': call
    }

# A8. Load existing nosz_clades.tsv for comparison
prior = {}
prior_file = DATA / 'nosz_clades.tsv'
if prior_file.exists():
    with open(prior_file) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            prior[row['gene_id']] = row

# A9. Write partA output
print('\nA9. Writing partA_clade_I_II.tsv...')
prior_to_cur = {}
for gene in nosz_genes:
    mag    = gene['mag']
    gid    = gene['gene_id']
    safe   = gid.replace('|', '_').replace('.', '_')
    pl     = placements.get(safe, {})
    mot    = motif_results.get(safe, {})
    pr     = prior.get(gid, {})
    prior_clade = pr.get('clade', 'N/A')
    tree_clade  = pl.get('clade', 'N/A')
    motif_call  = mot.get('motif_call', 'N/A')
    # Consensus: motif is deterministic; tree adds context
    if motif_call == 'L-NosZ_cladeIII':
        final_clade = 'III'
    elif tree_clade in ('I', 'II', 'III'):
        final_clade = tree_clade
    else:
        final_clade = prior_clade  # fall back to existing
    prior_to_cur[gid] = (prior_clade, final_clade)

fieldnames = ['gene_id', 'mag', 'prior_clade', 'final_clade',
              'tree_clade', 'nearest_ref', 'tree_dist', 'CuA', 'CuZ', 'motif_call']
with open(OUT / 'partA_clade_I_II.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    w.writeheader()
    for gene in nosz_genes:
        gid  = gene['gene_id']
        safe = gid.replace('|', '_').replace('.', '_')
        pl   = placements.get(safe, {})
        mot  = motif_results.get(safe, {})
        pr   = prior.get(gid, {})
        (prior_c, final_c) = prior_to_cur.get(gid, ('?', '?'))
        w.writerow({
            'gene_id':    gid,
            'mag':        gene['mag'],
            'prior_clade': prior_c,
            'final_clade': final_c,
            'tree_clade':  pl.get('clade', 'N/A'),
            'nearest_ref': pl.get('nearest_ref', 'N/A'),
            'tree_dist':   f"{pl.get('dist', float('nan')):.4f}" if pl else 'N/A',
            'CuA':         mot.get('CuA', 'N/A'),
            'CuZ':         mot.get('CuZ', 'N/A'),
            'motif_call':  mot.get('motif_call', 'N/A'),
        })

print(f'  Written: {OUT}/partA_clade_I_II.tsv')

# ═══════════════════════════════════════════════════════════════════════════════
# PART B — Scan all 276 bins for clade III L-NosZ
# ═══════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('PART B — Search all 276 bins for clade III L-NosZ')
print('=' * 70)

# B1. Concatenate all bin FASTAs with bin-prefixed contig names
print('\nB1. Building concatenated bin FASTA...')
all_bins_fna = OUT / 'all_bins.fna'
bin_fas = sorted(DREP.glob('*.fa'))
print(f'  Found {len(bin_fas)} bins in {DREP}')
if not all_bins_fna.exists():
    with open(all_bins_fna, 'w') as out_f:
        for fa in bin_fas:
            bin_name = fa.stem
            with open(fa) as in_f:
                for line in in_f:
                    if line.startswith('>'):
                        out_f.write(f'>{bin_name}::{line[1:]}')
                    else:
                        out_f.write(line)
    print(f'  Written: {all_bins_fna}')
else:
    print(f'  Using existing: {all_bins_fna}')

# B2. hmmsearch all bins with the He NosZ HMM (using --seg no for translated search)
# Since bins are nucleotide, we use --dnax option (DNA input, translated search)
print('\nB2. Running hmmsearch on all bins (translated 6-frame)...')
bins_tblout = OUT / 'partB_hmmsearch.tblout'
if not bins_tblout.exists():
    run([HMMSEARCH, '--cpu', '16', '--tblout', bins_tblout,
         '--noali', '--dnax',   # 6-frame DNA translation
         '-E', '1e-10',
         hmm, all_bins_fna])
else:
    print(f'  Using existing: {bins_tblout}')

# B3. Parse hmmsearch hits
print('\nB3. Parsing hmmsearch hits...')
bin_hits = {}  # hit_name -> score
with open(bins_tblout) as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.split()
        if len(p) < 10: continue
        name, score = p[0], float(p[5])
        bin_hits[name] = score

print(f'  Total bin hits: {len(bin_hits)}')

# B4. Extract hit sequences from all_bins_fna for protein extraction
# For nucleotide hits, extract the ORF and translate
# hmmsearch with --dnax reports the translated frame, but we need to get back the aa
# We'll use pyrodigal to call ORFs on hit contigs, then screen with HMM
print('\nB4. Extracting hit contigs for ORF calling...')
hit_contigs = set()
for name in bin_hits:
    # Hit name format from --dnax: contig/frame info
    # Take the base contig name
    hit_contigs.add(name.rsplit('/', 1)[0])

hit_fna = OUT / 'partB_hit_contigs.fna'
contig_seqs_all = {}
with open(all_bins_fna) as f:
    cur = None
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            cur = line[1:].split()[0]
            contig_seqs_all[cur] = []
        elif cur:
            contig_seqs_all[cur].append(line)

hit_seqs = {k: ''.join(v) for k, v in contig_seqs_all.items() if k in hit_contigs}
write_fasta(hit_fna, hit_seqs)
print(f'  {len(hit_seqs)} hit contigs extracted')

# B5. Call ORFs with pyrodigal on hit contigs, then screen with hmmsearch
print('\nB5. Calling ORFs on hit contigs...')
try:
    import pyrodigal
    HAS_PYRODIGAL = True
except ImportError:
    HAS_PYRODIGAL = False

hit_prots = {}
if HAS_PYRODIGAL and hit_seqs:
    orf_finder = pyrodigal.GeneFinder(meta=True)
    for contig_name, dna in hit_seqs.items():
        try:
            genes = orf_finder.find_genes(dna.encode())
            for i, gene in enumerate(genes):
                prot_name = f'{contig_name}_{i+1}'
                hit_prots[prot_name] = gene.translate()
        except Exception as e:
            print(f'  WARNING: pyrodigal failed on {contig_name}: {e}')
    print(f'  {len(hit_prots)} ORFs called from hit contigs')
else:
    if not HAS_PYRODIGAL:
        print('  pyrodigal not available; translating all 6 frames of hit contigs instead')
    # Fallback: translate all 6 frames, screen with HMM, take longest passing ORF
    for cname, dna in hit_seqs.items():
        for frame in range(3):
            for strand in ['+', '-']:
                aa = translate_dna(dna[frame:], strand)
                for j, frag in enumerate(aa.split('*')):
                    if len(frag) > 50:
                        prot_name = f'{cname}_f{frame}{strand}_orf{j}'
                        hit_prots[prot_name] = frag
    print(f'  {len(hit_prots)} 6-frame ORFs from hit contigs')

if hit_prots:
    hit_prots_faa = OUT / 'partB_hit_prots.faa'
    write_fasta(hit_prots_faa, hit_prots)

    # Confirm NosZ hits with HMM protein search
    print('\nB5b. Confirming NosZ hits with protein hmmsearch...')
    confirmed_tblout = OUT / 'partB_confirmed.tblout'
    run([HMMSEARCH, '--cpu', '8', '--tblout', confirmed_tblout,
         '--noali', '-E', '1e-5', hmm, hit_prots_faa])

    confirmed_hits = {}
    with open(confirmed_tblout) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.split()
            if len(p) < 10: continue
            confirmed_hits[p[0]] = float(p[5])
    print(f'  Confirmed NosZ ORFs: {len(confirmed_hits)}')

    # B6. Motif diagnosis on confirmed ORFs
    print('\nB6. Motif diagnosis on confirmed ORFs...')
    partB_motif_rows = []
    for prot_name, aa in hit_prots.items():
        if prot_name not in confirmed_hits:
            continue
        m = CUA.search(aa)
        tail = aa[m.end():] if m else aa
        c_hit = bool(CUZC.search(tail))
        l_hit = bool(CUZL.search(tail))
        if l_hit and not c_hit:
            call = 'L-NosZ_cladeIII'
        elif c_hit:
            call = 'C-NosZ_cladeI_II'
        else:
            call = 'ambiguous'
        bin_name = prot_name.split('::')[0] if '::' in prot_name else prot_name.rsplit('_', 1)[0]
        partB_motif_rows.append({
            'orf': prot_name, 'bin': bin_name,
            'CuA': 'CuA' if m else 'noCuA',
            'CuZ': 'GXHH' if (l_hit and not c_hit) else ('DXHH' if c_hit else '?'),
            'call': call, 'hmm_score': confirmed_hits[prot_name]
        })

    with open(OUT / 'partB_motif.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['orf','bin','CuA','CuZ','call','hmm_score'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(partB_motif_rows)
    print(f'  Written: {OUT}/partB_motif.tsv')

    # B7. L-NosZ candidates
    lnosz = [r for r in partB_motif_rows if r['call'] == 'L-NosZ_cladeIII']
    print(f'\n  L-NosZ (clade III) candidates: {len(lnosz)}')

    # B7b. Confirm L-NosZ: < 35% identity to Chee+Orellana C-NosZ via diamond
    lnosz_hits = []
    if lnosz:
        lnosz_faa = OUT / 'lnosz_candidates.faa'
        write_fasta(lnosz_faa, {r['orf']: hit_prots[r['orf']] for r in lnosz})
        chee_db = OUT / 'chee_diamond.dmnd'
        print('\nB7b. Building diamond db from Chee+Orellana C-NosZ...')
        run([DIAMOND, 'makedb', '--in', CHEE_FAA, '--db', chee_db, '--quiet'])
        diamond_out = OUT / 'lnosz_vs_cnosz.tsv'
        run([DIAMOND, 'blastp', '-q', lnosz_faa, '--db', chee_db,
             '--out', diamond_out, '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
             '--max-target-seqs', '1', '--quiet', '--evalue', '0.001'])
        cnosz_identity = {}
        with open(diamond_out) as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) >= 3:
                    cnosz_identity[p[0]] = float(p[2])
        for r in lnosz:
            pct = cnosz_identity.get(r['orf'], None)
            is_lnosz = pct is None or pct < 35.0
            lnosz_hits.append({**r,
                                'pct_id_to_CNosZ': f'{pct:.1f}%' if pct else 'no_hit',
                                'confirmed_LNosZ': 'YES' if is_lnosz else 'NO (>35% to C-NosZ)'})

        with open(OUT / 'partB_lnosz_hits.tsv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(lnosz_hits[0].keys()), delimiter='\t')
            w.writeheader()
            w.writerows(lnosz_hits)
        confirmed_lnosz = [r for r in lnosz_hits if r['confirmed_LNosZ'] == 'YES']
        print(f'  Confirmed L-NosZ (clade III, <35% to C-NosZ): {len(confirmed_lnosz)}')
else:
    print('  No hit contigs to process for Part B')
    partB_motif_rows = []
    lnosz = []
    lnosz_hits = []
    confirmed_lnosz = []

# ═══════════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('SUMMARY')
print('=' * 70)

changed = [(gid, pc, fc) for gid, (pc, fc) in prior_to_cur.items() if pc != fc and fc != '?']
print(f'\nPart A: {len(nosz_genes)} K00376 genes classified')
if changed:
    print(f'  Changes from prior call:')
    for gid, pc, fc in changed:
        print(f'    {gid}: {pc} → {fc}')
else:
    print('  No changes from prior clade I/II calls (or no placement available)')

motif_iii = [s for s, r in motif_results.items() if r['motif_call'] == 'L-NosZ_cladeIII']
print(f'\nMotif check (K00376 genes): {len(motif_iii)} with L-NosZ CuZ motif')

if hit_prots:
    confirmed_lnosz = locals().get('confirmed_lnosz', [])
    print(f'\nPart B: scanned {len(bin_fas)} bins')
    print(f'  hmmsearch hits: {len(bin_hits)}')
    print(f'  Confirmed NosZ ORFs: {len(confirmed_hits)}')
    print(f'  L-NosZ motif candidates: {len(lnosz)}')
    print(f'  Confirmed clade III (<35% to C-NosZ): {len(confirmed_lnosz)}')
else:
    print(f'\nPart B: scanned {len(bin_fas)} bins, {len(bin_hits)} raw hits')

summary_text = f"""# nosZ Clade Classification Summary

## qNOR / norZ status
- **K11188 (qNOR/norZ)**: completely absent from all 21 selected MAGs and from the
  entire 337-MAG dataset (no kofamscan or eggnog hits). These organisms do not carry
  quinol-dependent NOR.

## norB (K04561) status
- **K04561 (norB, cNOR large catalytic subunit)**: present in 10/21 selected MAGs
  at or above the KOfamScan threshold (120). All 10 norB-containing contigs were
  already in combined_ref.fa; norB has been added to the target gene BED and RPKM
  has been recomputed.
- norC (K02305) and norB (K04561) together constitute cNOR (cytochrome c-dependent
  NO reductase). The 4 MAGs with only norC in the BED previously undercounted cNOR
  activity; 6 additional MAGs now have cNOR coverage via norB only.

## Part A — K00376 nosZ clade I / II
- {len(nosz_genes)} K00376 genes; {len(changed)} changed from prior provisional call.
- {'No clade I/II changes.' if not changed else 'Changes: ' + ', '.join(f'{g} {p}→{f}' for g,p,f in changed)}
- CuZ motif check: {sum(1 for r in motif_results.values() if r['CuZ']=='DXHH')} DXHH (C-NosZ),
  {sum(1 for r in motif_results.values() if r['CuZ']=='GXHH')} GXHH (L-NosZ motif; investigate if present).

## Part B — Clade III L-NosZ in bins
- Scanned {len(bin_fas)} dereplicated bins with He 269NosZ HMM.
- Raw nucleotide hits: {len(bin_hits)}.
- Confirmed clade III L-NosZ (motif GXHH + <35% to C-NosZ): {'PRESENT — see partB_lnosz_hits.tsv' if locals().get('confirmed_lnosz') else 'NONE DETECTED'}.

## Figure integration
- The nosZ row in the heatmap uses clade calls from `out/partA_clade_I_II.tsv`.
- If clade calls changed, regenerate `gene_rpkm_per_sample_claded.tsv` and re-run
  `gene_ab_figure/04_plot_heatmap.py`.
"""
with open(OUT / 'summary.md', 'w') as f:
    f.write(summary_text)
print(f'\nWritten: {OUT}/summary.md')
print('\nAll done.')
