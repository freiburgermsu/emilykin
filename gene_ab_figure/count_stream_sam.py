#!/usr/bin/env python3
"""Stream minimap2 SAM on stdin; count distinct reads overlapping each target gene
and total distinct mapped reads (the RPKM 'per million mapped reads' denominator).

SAM variant of count_stream.py, so we can consume `minimap2 -ax map-ont` (base-level
alignment), matching the original 02_build_and_align.sh method — base-level alignment
resolves which homologous gene copy a read belongs to, unlike seed-chain-only PAF.

minimap2 emits a read's alignments consecutively, so distinct reads and per-gene
de-duplication are O(1) memory via the current query name. Secondary alignments are
skipped (--secondary=no leaves none); supplementary (chimeric) share the query name
and so are de-duplicated. Reference span is computed from the CIGAR.
"""
import sys, csv, json, argparse, re

ap = argparse.ArgumentParser()
ap.add_argument("--bed", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

genes = {}  # contig -> [(start, end, gene_id), ...]
with open(a.bed) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        genes.setdefault(row["prefixed_contig"], []).append(
            (int(row["start"]), int(row["end"]), row["gene_id"]))

_CIG = re.compile(r'(\d+)([MIDNSHP=X])')
_REF = set('MDN=X')          # CIGAR ops that consume the reference
def ref_len(cig):
    return 0 if cig == '*' else sum(int(n) for n, op in _CIG.findall(cig) if op in _REF)

counts = {}
total_mapped = 0
last_q = None
cur_genes = set()

def flush():
    for g in cur_genes:
        counts[g] = counts.get(g, 0) + 1

for line in sys.stdin:
    if line[0] == '@':
        continue
    p = line.split('\t')
    if len(p) < 6:
        continue
    flag = int(p[1])
    if flag & 0x4 or flag & 0x100:      # unmapped or secondary
        continue
    q = p[0]
    if q != last_q:
        if last_q is not None:
            flush(); total_mapped += 1
        cur_genes = set(); last_q = q
    glist = genes.get(p[2])
    if glist:
        ts = int(p[3]) - 1              # SAM POS is 1-based; BED is 0-based
        te = ts + ref_len(p[5])
        for gs, ge, gid in glist:
            if ts < ge and te > gs:
                cur_genes.add(gid)
if last_q is not None:
    flush(); total_mapped += 1

json.dump({"sample": a.sample, "total_mapped": total_mapped, "gene_counts": counts},
          open(a.out, "w"), indent=0)
sys.stderr.write(
    f"[count_sam] {a.sample}: mapped_reads={total_mapped:,} "
    f"target_gene_read_hits={sum(counts.values()):,} "
    f"genes_with_signal={sum(1 for v in counts.values() if v)}\n")
