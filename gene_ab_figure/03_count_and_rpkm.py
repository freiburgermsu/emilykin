#!/usr/bin/env python3
"""
Step 3: Count reads per target KO gene, calculate RPKM per sample and combined.

For merged KO groups (K00362/K00363, K02305/K00376), take the maximum RPKM
across all genes in the group for each MAG.

Output:
  data/gene_rpkm.tsv              - KO group x MAG (combined all samples)
  data/gene_rpkm_per_sample.tsv   - KO group x MAG per sample
  data/read_counts_raw.tsv        - raw counts (gene_id, sample, count)
"""
import csv, subprocess, os, sys, json
from collections import defaultdict

SAMTOOLS = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/88fdb48d4d745c55ec2cd90b407de422_/bin/samtools'
WORK     = '/scratch1/afreiburger/emilykin/gene_ab_figure'
DATA     = f'{WORK}/data'
SAMPLES  = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

KO_GROUPS = [
    ('nirB/D',    ['K00362', 'K00363']),
    ('nirS',      ['K15864']),
    ('nirK',      ['K00368']),
    ('norC/nosZ', ['K02305', 'K00376']),
]

def run_samtools(cmd, timeout=60):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.stdout.strip()
    except subprocess.TimeoutExpired:
        return ''

def count_reads(bam, contig, start, end):
    """Count reads overlapping [start,end) (0-based BED); samtools region is 1-based."""
    region = f"{contig}:{start+1}-{end}"
    out = run_samtools([SAMTOOLS, 'view', '-c', '-F', '4', bam, region])
    try:
        return int(out)
    except ValueError:
        return 0

def total_mapped(bam):
    """Get total mapped read count from flagstat."""
    out = run_samtools([SAMTOOLS, 'flagstat', bam], timeout=300)
    for line in out.splitlines():
        if 'mapped (' in line:
            return int(line.split()[0])
    return 1

# ── 1. Load target genes ─────────────────────────────────────────────────
print("Loading target genes...")
target_genes = []
with open(f'{DATA}/target_genes.bed') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        target_genes.append(row)

# Unique (contig, start, end) per gene_id
unique_genes = {}  # gene_id -> (contig, start, end, gene_len)
for row in target_genes:
    gid = row['gene_id']
    if gid not in unique_genes:
        s, e = int(row['start']), int(row['end'])
        unique_genes[gid] = (row['prefixed_contig'], s, e, max(e - s, 1))

# gene_id -> [(ko, mag)]
gene_info = defaultdict(list)
for row in target_genes:
    gene_info[row['gene_id']].append((row['ko'], row['mag']))

# Load MAG list
selected_mags = [l.strip() for l in open(f'{DATA}/selected_mag_list.txt') if l.strip()]

print(f"  {len(target_genes)} gene entries, {len(unique_genes)} unique genes")
print(f"  {len(selected_mags)} MAGs selected")

# ── 2. Check available BAMs ───────────────────────────────────────────────
per_sample_bams = {}
for s in SAMPLES:
    bam = f'{DATA}/aligned_{s}_sorted.bam'
    if os.path.exists(bam) and os.path.exists(bam + '.bai'):
        per_sample_bams[s] = bam
        print(f"  Found BAM: {s}")
    else:
        print(f"  WARNING: BAM not ready for {s}: {bam}")

if not per_sample_bams:
    print("ERROR: No BAM files found. Run alignment first.")
    sys.exit(1)

# ── 3. Get total mapped reads per sample ─────────────────────────────────
print("Getting total mapped reads per sample...")
per_sample_total = {}
for s, bam in per_sample_bams.items():
    n = total_mapped(bam)
    per_sample_total[s] = n
    print(f"  {s}: {n:,} mapped reads")

# ── 4. Count reads per gene per sample ───────────────────────────────────
print(f"Counting reads for {len(unique_genes)} genes × {len(per_sample_bams)} samples...")
gene_sample_counts = {}  # gene_id -> sample -> count

for idx, (gene_id, (contig, start, end, glen)) in enumerate(unique_genes.items()):
    gene_sample_counts[gene_id] = {}
    for s, bam in per_sample_bams.items():
        gene_sample_counts[gene_id][s] = count_reads(bam, contig, start, end)
    if (idx + 1) % 10 == 0:
        print(f"  {idx+1}/{len(unique_genes)} genes counted")

# ── 5. Calculate RPKM ────────────────────────────────────────────────────
def rpkm(count, glen, total):
    if total == 0 or glen == 0:
        return 0.0
    return (count / (glen / 1000.0)) / (total / 1e6)

# Gene-level RPKM per sample
gene_rpkm = {}  # gene_id -> sample -> rpkm
for gene_id, (contig, start, end, glen) in unique_genes.items():
    gene_rpkm[gene_id] = {}
    for s in per_sample_bams:
        cnt = gene_sample_counts[gene_id].get(s, 0)
        gene_rpkm[gene_id][s] = rpkm(cnt, glen, per_sample_total.get(s, 1))

# Combined RPKM (sum counts across samples, sum totals)
total_combined = sum(per_sample_total.values())
gene_rpkm_combined = {}
for gene_id, (contig, start, end, glen) in unique_genes.items():
    total_count = sum(gene_sample_counts[gene_id].get(s, 0) for s in per_sample_bams)
    gene_rpkm_combined[gene_id] = rpkm(total_count, glen, total_combined)

# ── 6. KO group × MAG matrix ─────────────────────────────────────────────
print("Building KO group × MAG RPKM matrices...")

# Combined
ko_mag_combined = defaultdict(lambda: defaultdict(float))
ko_mag_per_sample = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

for gene_id, pairs in gene_info.items():
    for ko, mag in pairs:
        for group_name, group_kos in KO_GROUPS:
            if ko in group_kos:
                # Max RPKM across genes in group for this MAG
                val_combined = gene_rpkm_combined.get(gene_id, 0)
                if val_combined > ko_mag_combined[group_name][mag]:
                    ko_mag_combined[group_name][mag] = val_combined
                for s in SAMPLES:
                    val_s = gene_rpkm.get(gene_id, {}).get(s, 0)
                    if val_s > ko_mag_per_sample[group_name][mag][s]:
                        ko_mag_per_sample[group_name][mag][s] = val_s

group_names = [g for g, _ in KO_GROUPS]

# ── 7. Write outputs ──────────────────────────────────────────────────────
# Combined RPKM matrix
with open(f'{DATA}/gene_rpkm.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko_group'] + sorted(selected_mags))
    for group in group_names:
        row = [group] + [f"{ko_mag_combined[group].get(m, 0):.4f}"
                         for m in sorted(selected_mags)]
        w.writerow(row)

# Per-sample RPKM
with open(f'{DATA}/gene_rpkm_per_sample.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko_group', 'mag'] + SAMPLES)
    for group in group_names:
        for mag in sorted(selected_mags):
            vals = [f"{ko_mag_per_sample[group][mag].get(s, 0):.4f}" for s in SAMPLES]
            w.writerow([group, mag] + vals)

# Raw counts
with open(f'{DATA}/read_counts_raw.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['gene_id', 'ko', 'mag', 'gene_len_bp'] + SAMPLES + ['combined_count'])
    for row in target_genes:
        gid = row['gene_id']
        glen = unique_genes.get(gid, ('','','',1))[3]
        counts_by_sample = [str(gene_sample_counts.get(gid, {}).get(s, 0)) for s in SAMPLES]
        combined = sum(gene_sample_counts.get(gid, {}).get(s, 0) for s in SAMPLES)
        w.writerow([gid, row['ko'], row['mag'], glen] + counts_by_sample + [combined])

print(f"\nWrote {DATA}/gene_rpkm.tsv")
print(f"Wrote {DATA}/gene_rpkm_per_sample.tsv")
print(f"Wrote {DATA}/read_counts_raw.tsv")

print("\n=== RPKM summary per KO group ===")
for group in group_names:
    vals = list(ko_mag_combined[group].values())
    n_nonzero = sum(1 for v in vals if v > 0)
    mx = max(vals) if vals else 0
    print(f"  {group:12s}: {n_nonzero:2d} MAGs with signal, max={mx:.2f} RPKM")

print(f"\nTotal mapped reads: {total_combined:,}")
print("Done. Next: run 04_plot_heatmap.py")
