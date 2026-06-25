#!/usr/bin/env python3
"""
Add ppk1 classification columns (ppk1_type, ppk1_species) to the root annotated
copy-number matrix ko_copy_number_matrix_nosZ.csv (nosZ clade cols already
reconciled to the Part-D tree call set). Idempotent. Verifies KO cells are
unchanged vs the raw matrix and reports nosZ + ppk1 coverage.
"""
import csv
from pathlib import Path

ROOT   = Path("/home/freiburger/Documents/EmilyKin")
MATRIX = ROOT / "ko_copy_number_matrix_nosZ.csv"
RAW    = ROOT / "ko_copy_number_matrix.csv"
PPK1   = ROOT / "ppk1_classify/ppk1_classification.tsv"
NEW    = ["ppk1_type", "ppk1_species"]

# ppk1 classification (the 5 Accumulibacter MAGs that carry ppk1)
ppk1 = {}
with PPK1.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ppk1[r["MAG"]] = (r.get("ppk1_type", ""), r.get("proposed_species", ""))

with MATRIX.open() as f:
    rows = list(csv.reader(f))
header, body = rows[0], rows[1:]
# drop any prior ppk1 columns (idempotent)
keep = [i for i, c in enumerate(header) if c not in NEW]
header = [header[i] for i in keep]
body   = [[r[i] for i in keep] for r in body]
mag_i  = header.index("MAG")

out_header = header + NEW
n_ppk1 = 0
with MATRIX.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(out_header)
    for r in body:
        t, sp = ppk1.get(r[mag_i], ("", ""))
        if t:
            n_ppk1 += 1
        w.writerow(r + [t, sp])

# ── verification ────────────────────────────────────────────────────────────
def load(p):
    with open(p) as fh:
        return list(csv.reader(fh))
raw = load(RAW); new = load(MATRIX)
new_h = new[0]
ko_cols = [c for c in raw[0] if c != "MAG"]
# reconstruct raw (MAG + KO cells) from the annotated file
idx = {c: new_h.index(c) for c in raw[0]}
recon = [raw[0]] + [[r[idx["MAG"]]] + [r[idx[c]] for c in ko_cols] for r in new[1:]]
print("KO cells identical to raw matrix:", raw == recon)
print(f"rows: {len(new)-1}  cols: {len(new_h)} (added {NEW})")
print(f"ppk1-classified MAGs populated: {n_ppk1}")
# report the ppk1 rows + nosZ tally
ci = new_h.index("nosZ_clade")
n1 = sum(1 for r in new[1:] if r[ci] == "I")
n2 = sum(1 for r in new[1:] if r[ci] == "II")
print(f"nosZ clade tally (should be 2 / 11): Clade I={n1}  Clade II={n2}")
tt, ts = new_h.index("ppk1_type"), new_h.index("ppk1_species")
print("\nppk1 rows in matrix:")
for r in new[1:]:
    if r[tt]:
        print(f"  {r[idx['MAG']]:14s} K00937={r[idx['K00937']]}  ppk1_type={r[tt]:22s} {r[ts]}")
