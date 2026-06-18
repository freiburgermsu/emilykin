#!/usr/bin/env python3
"""
Part C — search ALL assembly contigs (co-assembly + per-sample, INCLUDING the
unbinned fraction) for clade III L-NosZ. Higher-resolution extension of the
bins-only Part B (which only saw the 276 dereplicated MAGs).

Method is identical to Part B so results are directly comparable:
  pyrodigal meta ORFs  ->  He 269NosZ HMM (pyhmmer, E <= 1e-10)
  ->  CuA/CuZ motif (DXHH = C-NosZ clade I/II ; GXHH = L-NosZ clade III)
  ->  confirm each GXHH candidate by % identity to Chee+Orellana C-NosZ (<35% = L-NosZ)

Diamond is replaced by a Biopython local alignment for the (rare) confirm step,
since diamond is not installed locally; the <35% vs ~60-90% gap is robust to the
exact identity definition.

Runs entirely from the local uv venv (~/Documents/py_venv): pyrodigal + pyhmmer + Bio.
"""
import re, sys, time, csv, itertools
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import pyrodigal
import pyhmmer

HERE = Path('/home/freiburger/Documents/EmilyKin/clade_classify')
OUT  = HERE / 'contig_scan'; OUT.mkdir(exist_ok=True)
HMM_PATH = HERE / 'out' / 'He_269NosZ.hmm'
CHEE     = HERE / 'Chee_plus_Orellana_NosZ.prot.raw.faa'
META = Path('/home/freiburger/Documents/EmilyKin/meta/results')

SOURCES = [
    ('coasm', META / '08a_map_shortreads_co' / 'all_co_contigs.fasta'),
    ('CAN_1', META / '08_map_shortreads' / 'CAN_1' / 'CAN_1_contigs.fasta'),
    ('CAN_2', META / '08_map_shortreads' / 'CAN_2' / 'CAN_2_contigs.fasta'),
    ('CAN_3', META / '08_map_shortreads' / 'CAN_3' / 'CAN_3_contigs.fasta'),
    ('CAN_4', META / '08_map_shortreads' / 'CAN_4' / 'CAN_4_contigs.fasta'),
    ('CAN_5', META / '08_map_shortreads' / 'CAN_5' / 'CAN_5_contigs.fasta'),
]

EVALUE = 1e-10          # same inclusion threshold as Part B
CPUS   = 48
ORF_FAA = OUT / 'all_contigs_orf.faa'

CUA  = re.compile(r'C.{2}FC.{3}H.EM')   # CuA — all N2OR
CUZC = re.compile(r'D.HH')              # CuZ — C-NosZ (clade I/II)
CUZL = re.compile(r'G.HH')              # CuZ — L-NosZ (clade III)

GF = pyrodigal.GeneFinder(meta=True)    # meta mode = thread-safe, read-only profiles


def iter_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def chunked(it, n):
    it = iter(it)
    while True:
        c = list(itertools.islice(it, n))
        if not c:
            return
        yield c


def call_chunk(chunk):
    """ORF-call a list of (src, name, seq); return [(orf_id, aa), ...]."""
    res = []
    for src, name, seq in chunk:
        try:
            genes = GF.find_genes(seq.encode())
        except Exception:
            continue
        for i, gene in enumerate(genes):
            res.append((f'{src}::{name}_{i+1}', gene.translate()))
    return res


def main():
    t0 = time.time()
    for src, path in SOURCES:
        if not path.exists():
            print(f'FATAL: missing input {path}', flush=True); sys.exit(1)

    # ── Phase 1: ORF-call all contigs (pyrodigal meta), prefixed names ──────────
    print(f'[{time.time()-t0:.0f}s] Phase 1: ORF-calling all contigs (pyrodigal meta, {CPUS} threads)', flush=True)
    total = 0
    with open(ORF_FAA, 'w') as fo, ThreadPoolExecutor(max_workers=CPUS) as ex:
        for src, path in SOURCES:
            n_src = 0
            contigs = [(src, name, seq) for name, seq in iter_fasta(path)]
            futs = [ex.submit(call_chunk, ch) for ch in chunked(contigs, 500)]
            for fut in futs:
                for orf_id, aa in fut.result():
                    fo.write(f'>{orf_id}\n')
                    for j in range(0, len(aa), 80):
                        fo.write(aa[j:j+80] + '\n')
                    total += 1; n_src += 1
            del contigs
            print(f'  [{time.time()-t0:.0f}s] {src}: {n_src:,} ORFs (running total {total:,})', flush=True)
    print(f'[{time.time()-t0:.0f}s] Phase 1 done: {total:,} ORFs -> {ORF_FAA}', flush=True)

    # ── Phase 2: hmmsearch the He 269NosZ HMM over all ORFs ─────────────────────
    print(f'[{time.time()-t0:.0f}s] Phase 2: pyhmmer hmmsearch (He 269NosZ HMM, E<={EVALUE})', flush=True)
    aa_alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(str(HMM_PATH)) as hf:
        hmm = hf.read()
    with pyhmmer.easel.SequenceFile(str(ORF_FAA), digital=True, alphabet=aa_alpha) as sf:
        seqs = sf.read_block()
    print(f'  [{time.time()-t0:.0f}s] loaded {len(seqs):,} ORFs into memory; searching...', flush=True)

    hits = {}  # orf_id -> (score, evalue)
    for top in pyhmmer.hmmer.hmmsearch([hmm], seqs, cpus=CPUS, E=EVALUE, incE=EVALUE):
        for hit in top:
            ev = hit.evalue
            if ev <= EVALUE:
                nm = hit.name
                nm = nm.decode() if isinstance(nm, (bytes, bytearray)) else nm
                hits[nm] = (hit.score, ev)
    print(f'[{time.time()-t0:.0f}s] Phase 2 done: {len(hits):,} HMM hits at E<={EVALUE}', flush=True)

    # ── Phase 3: pull hit ORF sequences (single streaming pass) ─────────────────
    hit_aa = {}
    want = set(hits)
    cur, keep, buf = None, False, []
    with open(ORF_FAA) as f:
        for line in f:
            if line.startswith('>'):
                if keep:
                    hit_aa[cur] = ''.join(buf)
                cur = line[1:].split()[0]; keep = cur in want; buf = []
            elif keep:
                buf.append(line.strip())
    if keep:
        hit_aa[cur] = ''.join(buf)

    # ── Phase 4: CuA/CuZ motif diagnosis (same logic as partB_scan_bins.py) ─────
    rows = []
    for orf_id, (score, ev) in sorted(hits.items(), key=lambda x: -x[1][0]):
        aa = hit_aa.get(orf_id, '')
        src = orf_id.split('::', 1)[0]
        contig = orf_id.split('::', 1)[1].rsplit('_', 1)[0] if '::' in orf_id else orf_id
        m = CUA.search(aa)
        head = aa[:m.start()] if m else aa          # CuZ sits N-terminal of CuA
        c_hit, l_hit = bool(CUZC.search(head)), bool(CUZL.search(head))
        if l_hit and not c_hit:
            call = 'L-NosZ_cladeIII'
        elif c_hit:
            call = 'C-NosZ_cladeI_II'
        else:                                        # fallback: scan C-terminal of CuA
            tail = aa[m.end():] if m else aa
            if CUZC.search(tail):   call = 'C-NosZ_cladeI_II'
            elif CUZL.search(tail): call = 'L-NosZ_cladeIII'
            elif m:                 call = 'CuA_present_no_CuZ'
            else:                   call = 'no_CuA'
        cuz = 'GXHH' if call.startswith('L-NosZ') else ('DXHH' if 'cladeI_II' in call else '?')
        rows.append({'orf': orf_id, 'source': src, 'contig': contig,
                     'len_aa': len(aa), 'CuA': 'CuA' if m else 'noCuA', 'CuZ': cuz,
                     'call': call, 'hmm_score': round(score, 1), 'evalue': f'{ev:.1e}'})

    with open(OUT / 'partC_motif.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['orf','source','contig','len_aa','CuA','CuZ','call','hmm_score','evalue'], delimiter='\t')
        w.writeheader(); w.writerows(rows)

    lnosz = [r for r in rows if r['call'] == 'L-NosZ_cladeIII']
    print(f'[{time.time()-t0:.0f}s] Phase 4: motif done. GXHH (L-NosZ) candidates: {len(lnosz)}', flush=True)

    # ── Phase 5: confirm GXHH candidates by % identity to Chee+Orellana ─────────
    confirmed = []
    if lnosz:
        from Bio import SeqIO
        from Bio.Align import PairwiseAligner, substitution_matrices
        refs = list(SeqIO.parse(str(CHEE), 'fasta'))
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.mode = 'local'
        aligner.open_gap_score = -11; aligner.extend_gap_score = -1
        for r in lnosz:
            q = hit_aa[r['orf']]
            best_pid, best_ref = 0.0, None
            for ref in refs:
                s = str(ref.seq).replace('*', '')
                if not s:
                    continue
                aln = aligner.align(q, s)[0]
                # % identity over aligned columns
                ta, tb = aln.aligned
                ident = cols = 0
                qa, ra = aln.target, aln.query
                for (qs, qe), (rs, re_) in zip(ta, tb):
                    seg_q, seg_r = qa[qs:qe], ra[rs:re_]
                    cols += len(seg_q)
                    ident += sum(1 for a, b in zip(seg_q, seg_r) if a == b)
                pid = 100.0 * ident / cols if cols else 0.0
                if pid > best_pid:
                    best_pid, best_ref = pid, ref.id
            verdict = 'YES_L-NosZ' if best_pid < 35.0 else f'NO_C-NosZ(pid={best_pid:.1f}%)'
            confirmed.append({**r, 'best_chee_hit': best_ref or 'none',
                              'pct_id_to_CNosZ': f'{best_pid:.1f}%', 'confirmed_LNosZ': verdict})

    with open(OUT / 'partC_lnosz_hits.tsv', 'w', newline='') as f:
        cols = ['orf','source','contig','len_aa','CuA','CuZ','call','hmm_score','evalue',
                'best_chee_hit','pct_id_to_CNosZ','confirmed_LNosZ']
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for c in confirmed:
            w.writerow({k: c.get(k, '') for k in cols})

    # ── Summary ─────────────────────────────────────────────────────────────────
    n_dxhh = sum(1 for r in rows if r['CuZ'] == 'DXHH')
    n_gxhh = sum(1 for r in rows if r['CuZ'] == 'GXHH')
    n_other = len(rows) - n_dxhh - n_gxhh
    yes = [c for c in confirmed if c['confirmed_LNosZ'].startswith('YES')]
    print('\n===== PART C SUMMARY (all assembly contigs) =====', flush=True)
    print(f'Total ORFs searched:        {total:,}', flush=True)
    print(f'HMM hits (E<={EVALUE}):       {len(rows)}', flush=True)
    print(f'  DXHH (C-NosZ clade I/II):  {n_dxhh}', flush=True)
    print(f'  GXHH (L-NosZ candidate):   {n_gxhh}', flush=True)
    print(f'  no clear CuZ (fragments):  {n_other}', flush=True)
    print(f'Confirmed clade III L-NosZ (<35% id): {len(yes)}', flush=True)
    for c in (confirmed or []):
        print(f'    {c["orf"]}  CuZ={c["CuZ"]}  best={c["best_chee_hit"]}  id={c["pct_id_to_CNosZ"]}  -> {c["confirmed_LNosZ"]}', flush=True)
    if not yes:
        print('  => Clade III L-NosZ NOT DETECTED in any assembly contig (incl. unbinned).', flush=True)
    print(f'[{time.time()-t0:.0f}s] DONE. Outputs in {OUT}/', flush=True)


if __name__ == '__main__':
    main()
