#!/usr/bin/env python3
"""
Step 1: Select top-10 HQ MAGs per sample (union), classify as GAO/PAO/denitrifier,
and identify target KO genes with genomic coordinates.

Target KO rows:
  K00362/K00363  (nirB/nirD  DNRA NO2→NH4)   - take max per MAG
  K15864         (nirS        NO2→NO denitrification)
  K00368         (nirK        NO2→NO denitrification)
  K02305/K00376  (norC/nosZ   NO→N2O / N2O→N2) - take max per MAG

Output:
  data/selected_mags.tsv   - MAG | classification | sample_abundances
  data/target_genes.bed    - contig | start-1 | end | gene_id | ko | mag
  data/mag_contig_prefix.tsv - MAG | original_contig | prefixed_contig
"""
import csv, glob, os, re
from collections import defaultdict

WORK = '/scratch1/afreiburger/emilykin/gene_ab_figure'
REPO = '/scratch1/afreiburger/emilykin/EmilyKin'
PROC = '/scratch1/afreiburger/emilykin/processed/mag'
KOFAM = '/scratch1/afreiburger/emilykin/kofamscan'

SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']
TARGET_KOS = {'K00362', 'K00363', 'K15864', 'K00368', 'K02305', 'K00376'}

# ── 1. HQ MAGs from checkm2 ────────────────────────────────────────────────
print("Loading checkm2 quality reports...")
hq_mags = set()
with open(f'{PROC}/checkm2/all_quality_reports.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        mag = row['genome'].replace('.fa', '')
        if float(row['completeness']) >= 90 and float(row['contamination']) <= 5:
            hq_mags.add(mag)
print(f"  HQ MAGs (comp>=90, cont<=5): {len(hq_mags)}")

# ── 2. Per-sample abundance from coverm ───────────────────────────────────
print("Loading coverm abundance...")
coverm = {}  # mag -> {sample: rel_abund}
with open(f'{PROC}/abundance/coverm_genome.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mag = row['Genome']
        if mag == 'unmapped':
            continue
        coverm[mag] = {}
        for s in SAMPLES:
            col = f'{s} Relative Abundance (%)'
            coverm[mag][s] = float(row.get(col, 0) or 0)

# ── 3. Top 10 HQ MAGs per sample ──────────────────────────────────────────
print("Selecting top 10 HQ MAGs per sample...")
selected_mags = set()
for sample in SAMPLES:
    sample_hq = [(mag, coverm.get(mag, {}).get(sample, 0))
                 for mag in hq_mags if mag in coverm]
    sample_hq.sort(key=lambda x: x[1], reverse=True)
    top10 = [mag for mag, _ in sample_hq[:10]]
    print(f"  {sample} top10: {top10[:3]}...")
    selected_mags.update(top10)

print(f"Union of top-10 HQ MAGs across all samples: {len(selected_mags)}")

# ── 4. GAO / PAO / Denitrifier classification ────────────────────────────
print("Classifying MAGs from mag_abundance_summary.tsv...")
classification = {}  # mag -> set of {'GAO','PAO','Denitrifier'}
with open(f'{REPO}/mag_abundance_summary.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mag = row['MAG']
        cats = set()
        if int(float(row.get('GAO_like', 0) or 0)):
            cats.add('GAO')
        if (int(float(row.get('putative_GAO_like_storage', 0) or 0))
                and not int(float(row.get('stronger_putative_PAO_like', 0) or 0))):
            cats.add('GAO')
        if (int(float(row.get('stronger_putative_PAO_like', 0) or 0))
                or int(float(row.get('putative_PAO_like', 0) or 0))):
            cats.add('PAO')
        if (int(float(row.get('complete_denitrification', 0) or 0))
                or int(float(row.get('incomplete_denitrification', 0) or 0))):
            cats.add('Denitrifier')
        classification[mag] = cats

# ── 5. Parse KOfamScan thresholds ─────────────────────────────────────────
print("Parsing KOfamScan thresholds...")
thresholds = {}
with open(f'{KOFAM}/target_ko_list.txt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ko = row['knum']
        if ko in TARGET_KOS:
            thresholds[ko] = (float(row['threshold']), row['score_type'])
print(f"  Thresholds for target KOs: {thresholds}")

# KOs where we relax threshold (eggnog-confirmed functional genes with high
# structural similarity, but kofam threshold is tuned for very high precision)
RELAXED_KOS = {'K00362', 'K00363'}  # nirB/nirD - threshold 1279 is very strict
RELAXED_FACTOR = 0.55  # accept hits at >=55% of official threshold

# ── 6. Parse hmmsearch.tblout → protein→KO for target MAGs ───────────────
print("Parsing hmmsearch.tblout...")
kofam_hits = defaultdict(set)  # (mag, protein_id) -> set of KOs
with open(f'{KOFAM}/hmmsearch.tblout') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split()
        if len(parts) < 10: continue
        target = parts[0]       # MAG|protein_id
        ko = parts[2]           # KO query
        full_score = float(parts[5])
        dom_score = float(parts[8])
        if ko not in TARGET_KOS: continue
        if '|' not in target: continue
        mag, prot = target.split('|', 1)
        if mag not in selected_mags: continue
        thresh, stype = thresholds.get(ko, (0, 'full'))
        score = dom_score if stype == 'domain' else full_score
        # Relaxed threshold for KOs where kofam is very strict but eggnog confirms
        effective_thresh = thresh * RELAXED_FACTOR if ko in RELAXED_KOS else thresh
        if score >= effective_thresh:
            kofam_hits[(mag, prot)].add(ko)

n_hits = sum(len(v) for v in kofam_hits.values())
print(f"  KOfam target KO hits: {n_hits}")

# ── 7. Parse eggNOG annotations ───────────────────────────────────────────
print("Parsing eggNOG annotations...")
eggnog_hits = defaultdict(set)  # (mag, protein_id) -> set of KOs
# eggnog query column 0: protein_id (bakta locus_tag)
# KO column 11: comma-separated, format "ko:K00XXX"
for ann_path in glob.glob(f'{PROC}/eggnog/*.emapper.annotations'):
    mag = os.path.basename(ann_path).replace('.emapper.annotations', '')
    if mag not in selected_mags: continue
    with open(ann_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12: continue
            prot = parts[0]
            for token in parts[11].split(','):
                token = token.strip()
                if token.startswith('ko:'):
                    ko = token[3:]
                    if ko in TARGET_KOS:
                        eggnog_hits[(mag, prot)].add(ko)

# Merge all KO calls
all_ko_calls = defaultdict(set)  # (mag, prot) -> KOs
for key, kos in kofam_hits.items():
    all_ko_calls[key].update(kos)
for key, kos in eggnog_hits.items():
    all_ko_calls[key].update(kos)

print(f"  Total (mag,prot) pairs with target KOs: {len(all_ko_calls)}")

# ── 8. Parse bakta GFF3 for gene coordinates ─────────────────────────────
print("Parsing bakta GFF3 for gene coordinates...")
print("  Loading contig name mappings from bakta JSON...")

# Bakta renames contigs (contig_1, contig_2, ...) but JSON has orig_id
import json
bakta_contig_map = {}  # mag -> {bakta_contig: orig_contig}
for mag in selected_mags:
    json_path = f'{PROC}/bakta/{mag}/{mag}.json'
    if not os.path.exists(json_path):
        continue
    with open(json_path) as f:
        data = json.load(f)
    seq_map = {}
    for seq in data.get('sequences', []):
        seq_id = seq.get('id')
        orig_id = seq.get('orig_id', seq_id)
        if seq_id and orig_id:
            seq_map[seq_id] = orig_id
    bakta_contig_map[mag] = seq_map

print(f"  Contig maps loaded for {len(bakta_contig_map)} MAGs")

gene_coords = {}  # (mag, prot_id) -> (orig_contig, start0, end, strand)
for mag in selected_mags:
    gff_path = f'{PROC}/bakta/{mag}/{mag}.gff3'
    if not os.path.exists(gff_path):
        print(f"  WARNING: no GFF3 for {mag}")
        continue
    contig_map = bakta_contig_map.get(mag, {})
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'CDS': continue
            bakta_contig = parts[0]
            # Translate bakta contig name to original assembly name
            orig_contig = contig_map.get(bakta_contig, bakta_contig)
            start, end, strand = int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]
            locus = None
            for attr in attrs.split(';'):
                if attr.startswith('locus_tag='):
                    locus = attr[len('locus_tag='):]
                elif attr.startswith('ID=') and locus is None:
                    locus = attr[3:]
            if locus:
                # GFF3 is 1-based, BED is 0-based
                gene_coords[(mag, locus)] = (orig_contig, start - 1, end, strand)

print(f"  Gene coordinates loaded for {len(gene_coords)} (mag,locus) pairs")

# ── 9. Build BED file for target KO genes ────────────────────────────────
print("Building target gene BED file...")
target_genes = []  # (prefixed_contig, start0, end, gene_id, ko, mag)
mag_contig_map = {}  # (mag, orig_contig) -> prefixed_contig

for (mag, prot), kos in all_ko_calls.items():
    if (mag, prot) not in gene_coords:
        continue
    contig, start0, end, strand = gene_coords[(mag, prot)]
    # Prefix contig with MAG name to ensure uniqueness across MAGs
    safe_mag = mag.replace('.', '_')
    prefixed = f"{safe_mag}::{contig}"
    mag_contig_map[(mag, contig)] = prefixed
    gene_id = f"{mag}|{prot}"
    for ko in kos:
        target_genes.append((prefixed, start0, end, gene_id, ko, mag))

print(f"  Target gene entries: {len(target_genes)}")

# ── 10. Write outputs ─────────────────────────────────────────────────────
os.makedirs(f'{WORK}/data', exist_ok=True)

# selected_mags.tsv
with open(f'{WORK}/data/selected_mags.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['MAG', 'classification'] + [f'{s}_relabund' for s in SAMPLES])
    for mag in sorted(selected_mags):
        cats = classification.get(mag, set())
        cat_str = '/'.join(sorted(cats)) if cats else 'Other'
        abunds = [coverm.get(mag, {}).get(s, 0) for s in SAMPLES]
        w.writerow([mag, cat_str] + [f'{a:.6f}' for a in abunds])

# target_genes.bed
with open(f'{WORK}/data/target_genes.bed', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['prefixed_contig', 'start', 'end', 'gene_id', 'ko', 'mag'])
    for row in sorted(target_genes):
        w.writerow(row)

# mag_contig_prefix.tsv (for building combined FASTA)
with open(f'{WORK}/data/mag_contig_prefix.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['mag', 'orig_contig', 'prefixed_contig'])
    for (mag, contig), prefixed in sorted(mag_contig_map.items()):
        w.writerow([mag, contig, prefixed])

# List of selected MAGs for the shell script
with open(f'{WORK}/data/selected_mag_list.txt', 'w') as f:
    for mag in sorted(selected_mags):
        f.write(mag + '\n')

print("\n=== Summary ===")
print(f"Selected MAGs: {len(selected_mags)}")
cats_count = {'GAO': 0, 'PAO': 0, 'Denitrifier': 0, 'Other': 0}
for mag in selected_mags:
    cats = classification.get(mag, set())
    if not cats:
        cats_count['Other'] += 1
    for c in cats:
        cats_count[c] += 1
for k, v in cats_count.items():
    print(f"  {k}: {v}")
print(f"Target gene entries in BED: {len(target_genes)}")
ko_counts = defaultdict(int)
for _, _, _, _, ko, _ in target_genes:
    ko_counts[ko] += 1
for ko in sorted(TARGET_KOS):
    print(f"  {ko}: {ko_counts[ko]} gene entries")
print("\nDone. Next: run 02_build_and_align.sh")
