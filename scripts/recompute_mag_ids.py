#!/usr/bin/env python3
"""Recompute MAG iterativeIDs after the priority swap (a taxonomically AGREEING
mapped ASV now beats the GTDB-derived label; disagreeing ASVs are ignored).

ASV iterativeIDs are UNCHANGED, so this ONLY recomputes MAG ids; it updates
mag_iterativeID_old_to_new.{csv,json} (preserving the ORIGINAL 'old' columns from
the committed mapping) and scripts/_new_ids.json (mag_new/mag_case only). It does
NOT touch any ASV file or the ASV old->new mapping.

Run:  ~/Documents/py_venv/bin/python scripts/recompute_mag_ids.py
"""
import json, csv, re
from collections import defaultdict, Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HAND = Path(__file__).resolve().parent / "_new_ids.json"
HASDIGIT = re.compile(r"\d")
def numeric(s): return bool(HASDIGIT.search(str(s)))
def clean(v):
    s = str(v).strip(); return "" if s.lower() in ("", "nan", "none") else s
def norm(v):
    # Compensate for naming-convention differences (Ca_ prefixes, GTDB _A/_F/_m
    # suffixes, sensu_stricto, ...) by comparing the longest token of split('_').
    v = clean(v)
    return max(v.split("_"), key=len).lower() if v else ""
GTDB = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]
G2A = {"Domain": "Kingdom", "Phylum": "Phylum", "Class": "Class",
       "Order": "Order", "Family": "Family", "Genus": "Genus"}

taxonomy = json.loads((ROOT / "taxonomy.json").read_text())
ids      = json.loads((ROOT / "iterativeIDs.json").read_text())     # current ASV ids (unchanged)
hash2id  = {h: i for i, h in ids.items()}
relab    = json.loads((ROOT / "relative_abundance.json").read_text())
mean_ab  = {h: (sum(d.values()) / len(d) if d else 0.0) for h, d in relab.items()}
a2m      = json.loads((ROOT / "asv_to_mag_mapping.json").read_text())["asv_to_mag"]
mags     = json.loads((ROOT / "mag_abundance_summary.json").read_text())
NONWEAK  = {"species", "genus", "family"}
mapped   = defaultdict(list)
for h, rec in a2m.items():
    if rec.get("confidence_tier") in NONWEAK and rec.get("best_mag"):
        mapped[rec["best_mag"]].append(h)

def gtdb(m, r): return clean(mags.get(m, {}).get("taxonomy", {}).get(r, ""))
def deepest(m):
    ch = ("", "")
    for L in GTDB:
        v = gtdb(m, L)
        if v and not numeric(v): ch = (L, v)
    return ch

order = [r[0] for r in csv.reader(open(ROOT / "mag_abundance_by_day_intersection_with_ids.csv"))][1:]
order = [m for m in order if m]

mag_new, mag_case, mc = {}, {}, defaultdict(int)
for m in order:
    rr, rv = deepest(m)
    agree = []
    if rr and rr != "Domain" and rv:
        ar = G2A[rr]
        for h in mapped.get(m, []):
            av = clean(taxonomy.get(h, {}).get(ar, ""))
            if av and norm(av) == norm(rv):
                agree.append(h)
    if agree:                                              # most-abundant AGREEING ASV
        best = max(agree, key=lambda h: (mean_ab.get(h, 0.0), hash2id.get(h, "")))
        mag_new[m] = hash2id[best]; mag_case[m] = "asv_inherit"; continue
    g = gtdb(m, "Genus")
    if g and not numeric(g):
        mc[g] += 1; mag_new[m] = f"{g}.{mc[g]}_m"; mag_case[m] = "gtdb_genus"
    else:
        lbl = rv or "Bacteria"; mc[lbl] += 1
        mag_new[m] = f"{lbl}.{mc[lbl]}_m"; mag_case[m] = "gtdb_lowest"

# preserve ORIGINAL 'old' columns from the committed mapping
old_rows = {r["MAG"]: r for r in csv.DictReader(open(ROOT / "mag_iterativeID_old_to_new.csv"))}
with open(ROOT / "mag_iterativeID_old_to_new.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["MAG", "old_mag_iterativeID", "old_intersection_id", "new_mag_iterativeID", "assignment_case"])
    for m in order:
        o = old_rows.get(m, {})
        w.writerow([m, o.get("old_mag_iterativeID", ""), o.get("old_intersection_id", ""), mag_new[m], mag_case[m]])
(ROOT / "mag_iterativeID_old_to_new.json").write_text(json.dumps({m: mag_new[m] for m in order}))

H = json.loads(HAND.read_text())
H["mag_new"] = mag_new; H["mag_case"] = mag_case
HAND.write_text(json.dumps(H))

cc = Counter(mag_case.values())
dups = {v: c for v, c in Counter(mag_new.values()).items() if c > 1}
print("MAG cases:", dict(cc), "| total", len(mag_new))
print("duplicate MAG ids:", len(dups), dups)
print("numeric-root _m ids (should be 0):",
      sum(1 for m in order if mag_new[m].endswith("_m") and numeric(mag_new[m].rsplit('.',1)[0])))
inh = [m for m in order if mag_case[m] == "asv_inherit"]
print("inherited ids that are NOT valid current ASV ids:", sum(1 for m in inh if mag_new[m] not in ids))
cur = {r["MAG"]: r["iterativeID"] for r in csv.DictReader(open(ROOT / "mag_abundance_by_day_intersection_with_ids.csv"))}
flips = [m for m in order if cur[m] != mag_new[m]]
print(f"MAGs changed vs current committed: {len(flips)}")
for m in flips[:14]:
    print(f"   {m:16s} {cur[m]:22s} -> {mag_new[m]:22s} [{mag_case[m]}]")
