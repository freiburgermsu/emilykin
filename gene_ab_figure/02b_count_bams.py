#!/usr/bin/env python3
"""
Server-side count generator for 03c (21-MAG focal reference, NO re-alignment).

Produces data/counts_<sample>.json = {"total_mapped": N, "gene_counts": {gid: c}}
by counting each target_genes.bed interval over the EXISTING 21-MAG BAMs with the
same method that built the committed figure:

    samtools view -c -F 4 BAM contig:start+1-end     (records overlapping the gene)
    total_mapped = sum of idxstats 'mapped' over all refs
                 (verified == samtools flagstat 'mapped' == the committed constants)

This deliberately does NOT use the Mac's count_stream_sam.py / 02_align_local.sh,
which align to the 276-MAG community reference (handoff §0/§9: conserved genes lose
70-87% of reads on that ref). BAMs live in the working dir; outputs go to the repo
data dir so 03c (run from the repo) consumes them.
"""
import os, json, subprocess, csv

SAM   = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/88fdb48d4d745c55ec2cd90b407de422_/bin/samtools'
BAMDIR = '/scratch1/afreiburger/emilykin/gene_ab_figure/data'           # existing 21-MAG BAMs
REPO   = os.path.dirname(os.path.abspath(__file__))                      # repo gene_ab_figure
DATA   = os.path.join(REPO, 'data')
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']
BED    = os.path.join(DATA, 'target_genes.bed')

def idxstats_total(bam):
    out = subprocess.run([SAM, 'idxstats', bam], capture_output=True, text=True).stdout
    return sum(int(l.split('\t')[2]) for l in out.splitlines() if len(l.split('\t')) >= 3)

def count_region(bam, contig, start0, end):
    region = f'{contig}:{start0 + 1}-{end}'          # BED 0-based -> samtools 1-based
    out = subprocess.run([SAM, 'view', '-c', '-F', '4', bam, region],
                         capture_output=True, text=True).stdout.strip()
    return int(out) if out.isdigit() else 0

# unique genes (gene_id -> contig,start,end); a gene_id can have multiple KO rows
genes = {}
with open(BED) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        gid = r['gene_id']
        if gid not in genes:
            genes[gid] = (r['prefixed_contig'], int(r['start']), int(r['end']))
print(f'{len(genes)} unique genes in {BED}')

for s in SAMPLES:
    bam = os.path.join(BAMDIR, f'aligned_{s}_sorted.bam')
    if not os.path.exists(bam):
        raise SystemExit(f'missing BAM {bam}')
    total = idxstats_total(bam)
    gene_counts = {gid: count_region(bam, c, s0, e) for gid, (c, s0, e) in genes.items()}
    out = os.path.join(DATA, f'counts_{s}.json')
    json.dump({'sample': s, 'total_mapped': total, 'gene_counts': gene_counts},
              open(out, 'w'), indent=0)
    nz = sum(1 for v in gene_counts.values() if v)
    print(f'  {s}: total_mapped={total:,}  genes_with_signal={nz}/{len(genes)}  -> {out}')

print('Done.')
