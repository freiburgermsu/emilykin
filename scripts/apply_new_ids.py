#!/usr/bin/env python3
"""
apply_new_ids.py  —  Stage B: propagate the new iterativeIDs into every committed
data artifact that carries an ID, using the bijective maps from Stage A
(scripts/_new_ids.json).  Only ID columns/keys are rewritten; numeric data,
GTDB 'genus' helper columns and all other fields are preserved verbatim.

Run AFTER scripts/regen_iterative_ids.py.
"""
import json, csv, io
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
H = json.loads((Path(__file__).resolve().parent / "_new_ids.json").read_text())
asv_old2new = H["asv_old2new"]            # old ASV iid -> new ASV iid
mag_new     = H["mag_new"]                # MAG name -> new MAG id

def relabel_asv(iid):
    return asv_old2new.get(iid, iid)

report = []

# 1. abundances.csv  — header columns are ASV iterativeIDs ----------------------
p = ROOT / "abundances.csv"
with open(p) as fh:
    rd = csv.reader(fh); rows = list(rd)
hdr = rows[0]
new_hdr = [hdr[0]] + [relabel_asv(c) for c in hdr[1:]]
miss = sum(1 for c in hdr[1:] if c not in asv_old2new)
with open(p, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(new_hdr); w.writerows(rows[1:])
report.append(f"abundances.csv: relabeled {len(hdr)-1} columns ({miss} unmapped)")

# 2. mag_abundance_summary.tsv  — iterativeIDs (list) + mag_iterativeID ---------
p = ROOT / "mag_abundance_summary.tsv"
with open(p) as fh:
    rd = csv.reader(fh, delimiter="\t"); rows = list(rd)
hdr = rows[0]; idx = {c: i for i, c in enumerate(hdr)}
ii, mi, mc = idx["iterativeIDs"], idx["mag_iterativeID"], idx["MAG"]
n_mag = n_list = 0
for r in rows[1:]:
    mag = r[mc]
    if r[ii]:
        r[ii] = ";".join(relabel_asv(t) for t in r[ii].split(";") if t); n_list += 1
    if mag in mag_new:
        r[mi] = mag_new[mag]; n_mag += 1
with open(p, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerows(rows)
report.append(f"mag_abundance_summary.tsv: {n_mag} mag_iterativeID, {n_list} iterativeIDs lists")

# 3. mag_abundance_summary.json  — per-MAG iterativeIDs list + add mag_iterativeID
p = ROOT / "mag_abundance_summary.json"
mags = json.loads(p.read_text())
n = 0
for mag, rec in mags.items():
    if isinstance(rec.get("iterativeIDs"), list):
        rec["iterativeIDs"] = [relabel_asv(t) for t in rec["iterativeIDs"]]
    if mag in mag_new:
        rec["mag_iterativeID"] = mag_new[mag]; n += 1
p.write_text(json.dumps(mags))
report.append(f"mag_abundance_summary.json: {n} MAGs (iterativeIDs relabeled + mag_iterativeID set)")

# helper: rewrite a delimited table, setting one column from MAG-name lookup -----
def set_mag_col(path, delim, mag_col, id_col):
    p = ROOT / path
    with open(p) as fh:
        rd = csv.reader(fh, delimiter=delim); rows = list(rd)
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
        w = csv.writer(fh, delimiter=delim); w.writerows(rows)
    report.append(f"{path}: {n} {id_col} set ({miss} MAGs not in map)")

# 4-9 MAG-keyed tables ----------------------------------------------------------
set_mag_col("mag_abundance_by_day_intersection_with_ids.csv", ",", "MAG", "iterativeID")
set_mag_col("gene_ab_figure/data/taxonomy_labels.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("gene_ab_figure/offline_bundle/data/taxonomy_labels.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("gene_ab_figure/data/gene_copy_per_mag.tsv", "\t", "MAG", "iterativeID")
set_mag_col("ko_gene_abundance/ko_rpkm_combined.tsv", "\t", "MAG", "mag_iterativeID")
set_mag_col("ko_gene_abundance/ko_rpkm_per_sample.tsv", "\t", "mag", "iterativeID")
set_mag_col("ko_copy_number_matrix_nosZ.csv", ",", "MAG", "iterativeID")

# 10. gene_copy_per_mag.md  — markdown table, relabel iterativeID col by MAG -----
p = ROOT / "gene_ab_figure/data/gene_copy_per_mag.md"
lines = p.read_text().splitlines()
out = []
for ln in lines:
    if ln.startswith("|") and "---" not in ln:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) >= 2 and cells[0] in mag_new:
            cells[1] = mag_new[cells[0]]
            ln = "| " + " | ".join(cells) + " |"
    out.append(ln)
p.write_text("\n".join(out) + ("\n" if p.read_text().endswith("\n") else ""))
report.append("gene_copy_per_mag.md: iterativeID column relabeled")

# 11. nonzero_per_day.json  — {day: {asv_iid: value}} ---------------------------
p = ROOT / "nonzero_per_day.json"
d = json.loads(p.read_text())
d = {day: {relabel_asv(k): v for k, v in inner.items()} for day, inner in d.items()}
p.write_text(json.dumps(d))
report.append(f"nonzero_per_day.json: {len(d)} days relabeled")

# 12. significantly_connected_organisms.json  — list of asv iids ----------------
p = ROOT / "significantly_connected_organisms.json"
d = json.loads(p.read_text())
p.write_text(json.dumps([relabel_asv(x) for x in d]))
report.append(f"significantly_connected_organisms.json: {len(d)} relabeled")

# 13. FDR_passing_pairs.npy  — array of asv iids --------------------------------
p = ROOT / "FDR_passing_pairs.npy"
arr = np.load(p, allow_pickle=True)
np.save(p, np.array([relabel_asv(str(x)) for x in arr], dtype=object))
report.append(f"FDR_passing_pairs.npy: {len(arr)} relabeled")

print("\n".join(report))
