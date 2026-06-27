#!/usr/bin/env python3
"""Stage B (MAG-only): propagate the recomputed MAG iterativeIDs into the
MAG-keyed columns of every committed data file, keyed by MAG name -> mag_new.
ASV-keyed data (abundances, correlations, summary 'iterativeIDs' lists, nonzero,
FDR, network) is left untouched because ASV iterativeIDs did not change.
Run AFTER scripts/recompute_mag_ids.py."""
import json, csv
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
mag_new = json.loads((Path(__file__).resolve().parent / "_new_ids.json").read_text())["mag_new"]
report = []

def set_mag_col(path, delim, mag_col, id_col):
    p = ROOT / path
    with open(p) as fh:
        rows = list(csv.reader(fh, delimiter=delim))
    hdr = rows[0]; idx = {c: i for i, c in enumerate(hdr)}
    mc, ic = idx[mag_col], idx[id_col]
    n = miss = 0
    for r in rows[1:]:
        if not r:
            continue
        mag = r[mc]
        if mag in mag_new:
            r[ic] = mag_new[mag]; n += 1
        elif mag and mag != "TOTAL":
            miss += 1
    with open(p, "w", newline="") as fh:
        csv.writer(fh, delimiter=delim).writerows(rows)
    report.append(f"{path}: {n} {id_col} set ({miss} MAGs not in map)")

set_mag_col("mag_abundance_summary.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("mag_abundance_by_day_intersection_with_ids.csv", ",", "MAG", "iterativeID")
set_mag_col("gene_ab_figure/data/taxonomy_labels.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("gene_ab_figure/offline_bundle/data/taxonomy_labels.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("gene_ab_figure/data/gene_copy_per_mag.tsv", "\t", "MAG", "iterativeID")
set_mag_col("ko_gene_abundance/ko_rpkm_combined.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("ko_gene_abundance/ko_rpkm_per_sample.tsv", "\t", "mag", "iterativeID")
set_mag_col("ko_copy_number_matrix_nosZ.csv", ",", "MAG", "iterativeID")
set_mag_col("clade_classify/out/allbins_nosz_master.tsv", "\t", "MAG", "iterativeID")

# mag_abundance_summary.json: mag_iterativeID per MAG (leave 'iterativeIDs' ASV list)
p = ROOT / "mag_abundance_summary.json"
mags = json.loads(p.read_text()); n = 0
for mag, rec in mags.items():
    if mag in mag_new:
        rec["mag_iterativeID"] = mag_new[mag]; n += 1
p.write_text(json.dumps(mags))
report.append(f"mag_abundance_summary.json: {n} mag_iterativeID set")

print("\n".join(report))
