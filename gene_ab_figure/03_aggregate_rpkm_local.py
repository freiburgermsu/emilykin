#!/usr/bin/env python3
"""Step 3 (local): aggregate per-sample counts_<sample>.json into RPKM matrices.

RPKM(gene, sample) = (count / (gene_kb)) / (mapped_reads_million)
  - mapped_reads = distinct reads mapped to the 276-genome community reference
    for that sample (the 'per million mapped reads' denominator).
Combined RPKM sums counts across samples over summed mapped reads.
KO-group cell per MAG = max RPKM over all genes carrying any KO in the group
(handles the slashed KO rows K00362/K00363 and K02305/K00376).

Writes data/gene_rpkm.tsv, gene_rpkm_per_sample.tsv, read_counts_raw.tsv in the
exact format 04_plot_heatmap.py consumes.
"""
import csv, json, os
from collections import defaultdict

WORK = os.path.dirname(os.path.abspath(__file__))
DATA = f"{WORK}/data"
SAMPLES = ["CAN_1", "CAN_2", "CAN_3", "CAN_4", "CAN_5"]
KO_GROUPS = [
    ("nirB/D",    ["K00362", "K00363"]),
    ("nirS",      ["K15864"]),
    ("nirK",      ["K00368"]),
    ("norC/nosZ", ["K02305", "K00376"]),
]

# ── load target genes (coords + KO + MAG) ──────────────────────────────────
target_genes = []
unique_genes = {}                 # gene_id -> (contig, start, end, gene_len)
gene_info = defaultdict(list)     # gene_id -> [(ko, mag)]
with open(f"{DATA}/target_genes.bed") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        target_genes.append(row)
        gid = row["gene_id"]
        s, e = int(row["start"]), int(row["end"])
        unique_genes.setdefault(gid, (row["prefixed_contig"], s, e, max(e - s, 1)))
        gene_info[gid].append((row["ko"], row["mag"]))

selected_mags = [l.strip() for l in open(f"{DATA}/selected_mag_list.txt") if l.strip()]

# ── load per-sample counts ─────────────────────────────────────────────────
per_total = {}                          # sample -> mapped reads
gene_counts = defaultdict(dict)         # gene_id -> sample -> count
for s in SAMPLES:
    p = f"{DATA}/counts_{s}.json"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} — alignment for {s} not finished")
    d = json.load(open(p))
    per_total[s] = d["total_mapped"]
    for gid, c in d["gene_counts"].items():
        gene_counts[gid][s] = c
    print(f"  {s}: {d['total_mapped']:,} mapped reads, "
          f"{sum(d['gene_counts'].values()):,} target-gene read hits")

def rpkm(count, glen, total):
    if not total or not glen:
        return 0.0
    return (count / (glen / 1000.0)) / (total / 1e6)

# ── gene-level RPKM ─────────────────────────────────────────────────────────
gene_rpkm = {}            # gene_id -> sample -> rpkm
gene_rpkm_combined = {}   # gene_id -> rpkm
total_combined = sum(per_total.values())
for gid, (c, s0, e0, glen) in unique_genes.items():
    gene_rpkm[gid] = {s: rpkm(gene_counts[gid].get(s, 0), glen, per_total[s]) for s in SAMPLES}
    tot = sum(gene_counts[gid].get(s, 0) for s in SAMPLES)
    gene_rpkm_combined[gid] = rpkm(tot, glen, total_combined)

# ── KO-group × MAG (max over group genes) ──────────────────────────────────
ko_mag_combined = defaultdict(lambda: defaultdict(float))
ko_mag_per_sample = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
for gid, pairs in gene_info.items():
    for ko, mag in pairs:
        for gname, gkos in KO_GROUPS:
            if ko in gkos:
                if gene_rpkm_combined[gid] > ko_mag_combined[gname][mag]:
                    ko_mag_combined[gname][mag] = gene_rpkm_combined[gid]
                for s in SAMPLES:
                    if gene_rpkm[gid][s] > ko_mag_per_sample[gname][mag][s]:
                        ko_mag_per_sample[gname][mag][s] = gene_rpkm[gid][s]

group_names = [g for g, _ in KO_GROUPS]

with open(f"{DATA}/gene_rpkm.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ko_group"] + sorted(selected_mags))
    for g in group_names:
        w.writerow([g] + [f"{ko_mag_combined[g].get(m, 0):.4f}" for m in sorted(selected_mags)])

with open(f"{DATA}/gene_rpkm_per_sample.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ko_group", "mag"] + SAMPLES)
    for g in group_names:
        for m in sorted(selected_mags):
            w.writerow([g, m] + [f"{ko_mag_per_sample[g][m].get(s, 0):.4f}" for s in SAMPLES])

with open(f"{DATA}/read_counts_raw.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id", "ko", "mag", "gene_len_bp"] + SAMPLES + ["combined_count"])
    for row in target_genes:
        gid = row["gene_id"]
        glen = unique_genes[gid][3]
        cs = [str(gene_counts[gid].get(s, 0)) for s in SAMPLES]
        w.writerow([gid, row["ko"], row["mag"], glen] + cs +
                   [sum(gene_counts[gid].get(s, 0) for s in SAMPLES)])

print(f"\nTotal mapped reads (all samples): {total_combined:,}")
print("=== RPKM per KO group (combined) ===")
for g in group_names:
    vals = list(ko_mag_combined[g].values())
    print(f"  {g:11s}: {sum(1 for v in vals if v>0):2d} MAGs with signal, "
          f"max={max(vals) if vals else 0:.2f} RPKM")
print("\nWrote gene_rpkm.tsv, gene_rpkm_per_sample.tsv, read_counts_raw.tsv")
