#!/usr/bin/env python3
"""Stream minimap2 PAF on stdin; count reads overlapping each target gene and
total distinct mapped reads (the RPKM 'per million mapped reads' denominator).

minimap2 emits all alignments of a read consecutively, so distinct reads and
per-read/per-gene de-duplication are done with O(1) memory by tracking the
current query name. PAF target coords (cols 8/9) are 0-based half-open, matching
the BED gene intervals (start = bakta_start-1, end = bakta_end).
"""
import sys, csv, json, argparse

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

counts = {}            # gene_id -> reads overlapping
total_mapped = 0       # distinct mapped reads
last_q = None
cur_genes = set()

def flush():
    for g in cur_genes:
        counts[g] = counts.get(g, 0) + 1

for line in sys.stdin:
    p = line.split("\t")
    if len(p) < 9:
        continue
    q = p[0]
    if q != last_q:
        if last_q is not None:
            flush()
        cur_genes = set()
        total_mapped += 1
        last_q = q
    glist = genes.get(p[5])
    if glist:
        ts, te = int(p[7]), int(p[8])
        for gs, ge, gid in glist:
            if ts < ge and te > gs:
                cur_genes.add(gid)
if last_q is not None:
    flush()

json.dump({"sample": a.sample, "total_mapped": total_mapped,
           "gene_counts": counts}, open(a.out, "w"), indent=0)
sys.stderr.write(
    f"[count_stream] {a.sample}: mapped_reads={total_mapped:,} "
    f"target_gene_read_hits={sum(counts.values()):,} "
    f"genes_with_signal={sum(1 for v in counts.values() if v)}\n")
