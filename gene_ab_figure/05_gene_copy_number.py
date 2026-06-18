#!/usr/bin/env python3
"""
Step 5: Compute gene copy number per MAG per sample.

Copy number = (reads_on_gene / gene_len_bp) / (reads_on_genome / genome_len_bp)

This normalizes each gene's coverage by its MAG's overall genome coverage,
yielding how many copies of the gene exist per genome copy in the community.

For merged KO groups (same as 03_count_and_rpkm.py), take the maximum copy
number across all genes in the group for each MAG × sample.

Outputs:
  data/gene_copy_number.tsv              - KO group × MAG (combined across samples)
  data/gene_copy_number_per_sample.tsv   - KO group × MAG × sample
  data/mag_genome_coverage.tsv           - genome coverage stats per MAG per sample
"""
import csv, subprocess, sys
from collections import defaultdict

SAMTOOLS = '/scratch1/afreiburger/emilykin/processed/.snakemake_envs/88fdb48d4d745c55ec2cd90b407de422_/bin/samtools'
DATA     = '/scratch1/afreiburger/emilykin/gene_ab_figure/data'
SAMPLES  = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

KO_GROUPS = [
    ('nirB/D',     ['K00362', 'K00363']),
    ('nirS',       ['K15864']),
    ('nirK',       ['K00368']),
    ('norBC/nosZ', ['K02305', 'K04561', 'K00376']),
]

# ── 1. Run samtools idxstats on each BAM ─────────────────────────────────
# idxstats columns: contig_name, contig_len, mapped_reads, unmapped_reads

print("Running samtools idxstats on BAMs...")

# contig_name -> {'len': int, sample -> int}
contig_stats = defaultdict(lambda: {'len': 0})

for sample in SAMPLES:
    bam = f'{DATA}/aligned_{sample}_sorted.bam'
    result = subprocess.run(
        [SAMTOOLS, 'idxstats', bam],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        print(f"  ERROR running idxstats on {bam}: {result.stderr.strip()}")
        sys.exit(1)
    n_contigs = 0
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        contig_name = parts[0]
        if contig_name == '*':  # unmapped bucket
            continue
        contig_len   = int(parts[1])
        mapped_reads = int(parts[2])
        contig_stats[contig_name]['len'] = contig_len
        contig_stats[contig_name][sample] = mapped_reads
        n_contigs += 1
    print(f"  {sample}: {n_contigs} contigs indexed")

# ── 2. Aggregate to MAG-level genome stats ────────────────────────────────
# MAG name is the part before "::" in the prefixed contig name

print("Aggregating contig stats to MAG level...")

mag_genome_len    = defaultdict(int)   # mag -> total bp
mag_mapped_reads  = defaultdict(lambda: defaultdict(int))  # mag -> sample -> reads

for contig_name, stats in contig_stats.items():
    if '::' in contig_name:
        mag = contig_name.split('::')[0].replace('_bin_', '.bin.').replace('_', '.', 1)
        # Reconstruct MAG ID: "CAN_1_bin_210" -> "CAN_1_bin.210"
        # Format from BAM header: CAN_1_bin_210, coasm_bin_185, etc.
        # Need to match the MAG IDs used in read_counts_raw.tsv: CAN_1_bin.210, coasm_bin.185
        prefix = contig_name.split('::')[0]  # e.g. CAN_1_bin_210
        # Find last underscore before the bin number and replace with "."
        # Pattern: everything up to the last "_" in "bin_NNN"
        last_underscore = prefix.rfind('_')
        mag = prefix[:last_underscore] + '.' + prefix[last_underscore+1:]
        # e.g. CAN_1_bin.210, coasm_bin.185
    else:
        mag = contig_name
    mag_genome_len[mag] += stats['len']
    for sample in SAMPLES:
        mag_mapped_reads[mag][sample] += stats.get(sample, 0)

print(f"  {len(mag_genome_len)} MAGs found in reference")
for mag in sorted(mag_genome_len):
    print(f"    {mag}: {mag_genome_len[mag]:,} bp  "
          f"reads={[mag_mapped_reads[mag][s] for s in SAMPLES]}")

# ── 3. Load gene-level raw counts ────────────────────────────────────────
print("\nLoading gene read counts...")
genes = []
with open(f'{DATA}/read_counts_raw.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        genes.append({
            'gene_id': row['gene_id'],
            'ko':      row['ko'],
            'mag':     row['mag'],
            'gene_len': int(row['gene_len_bp']),
            **{s: int(row[s]) for s in SAMPLES},
        })
print(f"  {len(genes)} gene entries loaded")

# ── 4. Compute per-gene copy numbers ─────────────────────────────────────
print("Computing gene copy numbers...")

def copy_number(reads_gene, gene_len, reads_genome, genome_len):
    if reads_genome == 0 or genome_len == 0 or gene_len == 0:
        return 0.0
    gene_depth   = reads_gene   / gene_len
    genome_depth = reads_genome / genome_len
    return gene_depth / genome_depth

gene_cn = []   # list of dicts: gene_id, ko, mag, sample -> cn, combined_cn
for g in genes:
    mag         = g['mag']
    gene_len    = g['gene_len']
    genome_len  = mag_genome_len.get(mag, 0)

    per_sample_cn = {}
    for s in SAMPLES:
        reads_gene   = g[s]
        reads_genome = mag_mapped_reads[mag][s]
        per_sample_cn[s] = copy_number(reads_gene, gene_len, reads_genome, genome_len)

    # Combined: use summed reads across all samples
    total_gene_reads   = sum(g[s] for s in SAMPLES)
    total_genome_reads = sum(mag_mapped_reads[mag][s] for s in SAMPLES)
    combined_cn = copy_number(total_gene_reads, gene_len, total_genome_reads, genome_len)

    gene_cn.append({
        'gene_id':    g['gene_id'],
        'ko':         g['ko'],
        'mag':        mag,
        'gene_len':   gene_len,
        'genome_len': genome_len,
        'combined_cn': combined_cn,
        **{f'{s}_cn': per_sample_cn[s] for s in SAMPLES},
    })

# ── 5. Aggregate to KO group × MAG (max across genes in group) ────────────
print("Aggregating to KO group × MAG matrices...")

all_mags = sorted(set(g['mag'] for g in gene_cn))

ko_mag_combined = defaultdict(lambda: defaultdict(float))
ko_mag_per_sample = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

for g in gene_cn:
    ko  = g['ko']
    mag = g['mag']
    for group_name, group_kos in KO_GROUPS:
        if ko in group_kos:
            if g['combined_cn'] > ko_mag_combined[group_name][mag]:
                ko_mag_combined[group_name][mag] = g['combined_cn']
            for s in SAMPLES:
                val = g[f'{s}_cn']
                if val > ko_mag_per_sample[group_name][mag][s]:
                    ko_mag_per_sample[group_name][mag][s] = val

group_names = [g for g, _ in KO_GROUPS]

# ── 6. Write outputs ──────────────────────────────────────────────────────

# Combined copy number matrix
with open(f'{DATA}/gene_copy_number.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko_group'] + all_mags)
    for group in group_names:
        row = [group] + [f"{ko_mag_combined[group].get(m, 0):.6f}" for m in all_mags]
        w.writerow(row)
print(f"Wrote {DATA}/gene_copy_number.tsv")

# Per-sample copy number (long format: ko_group, mag, CAN_1..CAN_5)
with open(f'{DATA}/gene_copy_number_per_sample.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['ko_group', 'mag'] + SAMPLES)
    for group in group_names:
        for mag in all_mags:
            vals = [f"{ko_mag_per_sample[group][mag].get(s, 0):.6f}" for s in SAMPLES]
            w.writerow([group, mag] + vals)
print(f"Wrote {DATA}/gene_copy_number_per_sample.tsv")

# MAG genome coverage table
with open(f'{DATA}/mag_genome_coverage.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['mag', 'genome_len_bp'] + [f'{s}_mapped_reads' for s in SAMPLES]
               + [f'{s}_avg_depth' for s in SAMPLES])
    for mag in sorted(mag_genome_len):
        gl = mag_genome_len[mag]
        reads_per_sample = [mag_mapped_reads[mag][s] for s in SAMPLES]
        # avg depth ≈ reads * avg_read_len / genome_len  — use reads/bp as proxy (no read len)
        # We'll just output reads/genome_len as a coverage proxy
        depth_per_sample = [f"{mag_mapped_reads[mag][s]/gl:.4f}" for s in SAMPLES]
        w.writerow([mag, gl] + reads_per_sample + depth_per_sample)
print(f"Wrote {DATA}/mag_genome_coverage.tsv")

# ── 7. Summary ────────────────────────────────────────────────────────────
print("\n=== Copy number summary per KO group ===")
for group in group_names:
    vals = [ko_mag_combined[group].get(m, 0) for m in all_mags]
    nonzero = [(m, v) for m, v in zip(all_mags, vals) if v > 0]
    if nonzero:
        mx_mag, mx_val = max(nonzero, key=lambda x: x[1])
        print(f"  {group:15s}: {len(nonzero):2d} MAGs with signal, "
              f"max={mx_val:.4f} (in {mx_mag})")
    else:
        print(f"  {group:15s}: no signal")

print("\n=== Full combined copy number matrix ===")
header = ['ko_group'] + all_mags
print('\t'.join(header))
for group in group_names:
    vals = [f"{ko_mag_combined[group].get(m, 0):.4f}" for m in all_mags]
    print('\t'.join([group] + vals))

print("\nDone.")
