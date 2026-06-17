#!/usr/bin/env python3
"""Step 1 of Part B: call ORFs on all 276 bins with pyrodigal (meta mode)."""
import csv, sys
from pathlib import Path
try:
    import pyrodigal
except ImportError:
    sys.exit('ERROR: pyrodigal not available in this env. Run with env1 python.')

DREP = Path('/scratch1/afreiburger/emilykin/processed/mag/drep/dereplicated_genomes')
OUT  = Path('/scratch1/afreiburger/emilykin/EmilyKin/clade_classify/out')
OUT.mkdir(exist_ok=True)
OUT_FAA = OUT / 'all_bins_orf.faa'

bin_fas = sorted(DREP.glob('*.fa'))
print(f'Calling ORFs on {len(bin_fas)} bins with pyrodigal (meta mode)...')
finder = pyrodigal.GeneFinder(meta=True)

total = 0
with open(OUT_FAA, 'w') as out_f:
    for fa in bin_fas:
        bn = fa.stem
        seqs = {}; cur = None
        with open(fa) as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    cur = line[1:].split()[0]; seqs[cur] = []
                elif cur:
                    seqs[cur].append(line)
        seqs = {k: ''.join(v).encode() for k, v in seqs.items()}
        for contig, dna in seqs.items():
            try:
                genes = finder.find_genes(dna)
                for i, gene in enumerate(genes):
                    prot_id = f'{bn}::{contig}_{i+1}'
                    aa = gene.translate()
                    out_f.write(f'>{prot_id}\n')
                    for j in range(0, len(aa), 80):
                        out_f.write(aa[j:j+80] + '\n')
                    total += 1
            except Exception as e:
                pass  # skip problematic contigs
        if (bin_fas.index(fa) + 1) % 50 == 0:
            print(f'  {bin_fas.index(fa)+1}/{len(bin_fas)} bins processed, {total:,} ORFs so far')

print(f'Done: {total:,} ORFs written to {OUT_FAA}')
