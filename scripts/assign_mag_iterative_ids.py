#!/usr/bin/env python3
"""
assign_mag_iterative_ids.py

Reassign the ``iterativeID`` column of
``mag_abundance_by_day_intersection_with_ids.csv`` (one iterativeID per MAG)
using the current nomenclature.

Assignment rules (in priority order):

  case 1 — GTDB genus root                       =>  ``{Genus}.{N}_m``
      The GTDB Genus is defined and non-numeric (contains no digit).

  case 2 — inherit a mapped ASV's iterativeID     =>  ``{ASV_iterativeID}``  (no suffix)
      Reached when case 1 does not apply.  Among the ASVs mapped to this MAG at
      species/genus/family confidence (the ``weak`` tier is excluded), keep those
      whose taxonomy agrees with the MAG's GTDB taxonomy at the DEEPEST GTDB rank
      that is defined and non-numeric (Phylum..Genus; Domain is never a valid
      match rank).  Of the survivors, inherit the iterativeID of the one with the
      highest mean relative abundance.  The MAG and that representative ASV then
      share a label (e.g. ``Bdellovibrio.13``).

  case 3 — lowest GTDB label root                 =>  ``{label}.{N}_m``
      Reached when GTDB does not resolve the genus (undefined/numeric) AND no
      mapped ASV agrees.  ``label`` = the deepest GTDB rank value (Domain..Genus)
      that is defined and non-numeric.

  ``_m`` is appended AFTER the per-root running count; case-2 (inherited) IDs
  carry no suffix.  Per-root counters are shared across cases 1 and 3.

GTDB taxonomy is read from the authoritative ``mag_abundance_summary.json``.
ASV iterativeIDs / taxonomy come from the (regenerated) ``iterativeIDs.json`` and
``taxonomy.json``.  Only the ``iterativeID`` column of the CSV is rewritten.

Run:  ~/Documents/py_venv/bin/python scripts/assign_mag_iterative_ids.py
"""
from collections import defaultdict, Counter
from pathlib import Path
import json, re

ROOT = Path(__file__).resolve().parent.parent

IDS_JSON  = ROOT / "iterativeIDs.json"            # {iterativeID: asv_hash}  (current scheme)
TAXONOMY  = ROOT / "taxonomy.json"                # {asv_hash: {rank: value}}
ASV_MAG   = ROOT / "asv_to_mag_mapping.json"      # {asv_hash: {best_mag, confidence_tier, ...}}
MAG_SUM   = ROOT / "mag_abundance_summary.json"   # {mag: {taxonomy: {...}, ...}}
RELAB     = ROOT / "relative_abundance.json"      # {asv_hash: {sample: rel_abundance}}
CSV       = ROOT / "mag_abundance_by_day_intersection_with_ids.csv"

NONWEAK = {"species", "genus", "family"}
GTDB_LEVELS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]   # Species excluded
GTDB2ASV_RANK = {"Domain": "Kingdom", "Phylum": "Phylum", "Class": "Class",
                 "Order": "Order", "Family": "Family", "Genus": "Genus"}
HASDIGIT = re.compile(r"\d")

def numeric(s):
    return bool(HASDIGIT.search(str(s)))
def _clean(v):
    s = str(v).strip()
    return "" if s.lower() in ("", "nan", "none") else s

ids       = json.loads(IDS_JSON.read_text())          # iid -> hash
hash2id   = {h: i for i, h in ids.items()}
taxonomy  = json.loads(TAXONOMY.read_text())          # hash -> {rank:val}
a2m       = json.loads(ASV_MAG.read_text())["asv_to_mag"]
mags      = json.loads(MAG_SUM.read_text())
relab     = json.loads(RELAB.read_text())
mean_ab   = {h: (sum(d.values()) / len(d) if d else 0.0) for h, d in relab.items()}

mag_to_hashes = defaultdict(list)
for h, rec in a2m.items():
    if rec.get("confidence_tier") in NONWEAK and rec.get("best_mag"):
        mag_to_hashes[rec["best_mag"]].append(h)

def gtdb(mag, rank):
    return _clean(mags.get(mag, {}).get("taxonomy", {}).get(rank, ""))

def deepest_nonnum_gtdb(mag):
    chosen = ("", "")
    for L in GTDB_LEVELS:
        v = gtdb(mag, L)
        if v and not numeric(v):
            chosen = (L, v)
    return chosen

def assign(mag, counter):
    g = gtdb(mag, "Genus")
    if g and not numeric(g):                                    # case 1
        counter[g] += 1
        return f"{g}.{counter[g]}_m", "1_gtdb_genus"
    r_rank, r_val = deepest_nonnum_gtdb(mag)
    agree = []
    if r_rank and r_rank != "Domain" and r_val:                 # case 2
        ar = GTDB2ASV_RANK[r_rank]
        for h in mag_to_hashes.get(mag, []):
            av = _clean(taxonomy.get(h, {}).get(ar, ""))
            if av and av.lower() == r_val.lower():
                agree.append(h)
    if agree:
        best = max(agree, key=lambda h: (mean_ab.get(h, 0.0), hash2id.get(h, "")))
        return hash2id[best], "2_asv_inherit"
    lbl = r_val if r_val else "Bacteria"                        # case 3
    counter[lbl] += 1
    return f"{lbl}.{counter[lbl]}_m", "3_lowest_gtdb"

# ── compute in CSV row order ──────────────────────────────────────────────────
lines = CSV.read_text().splitlines(keepends=True)
header, body = lines[0], lines[1:]
assigned, kinds, counter = {}, {}, defaultdict(int)
for line in body:
    mag = line.split(",", 1)[0]
    if not mag:
        continue
    iid, kind = assign(mag, counter)
    assigned[mag], kinds[mag] = iid, kind

out = [header]
for line in body:
    nl = "\n" if line.endswith("\n") else ""
    core = line[:-len(nl)] if nl else line
    if not core:
        out.append(line); continue
    mag, _old, rest = core.split(",", 2)
    out.append(f"{mag},{assigned[mag]},{rest}{nl}")
CSV.write_text("".join(out))

cc = Counter(kinds.values())
dups = {v: c for v, c in Counter(assigned.values()).items() if c > 1}
print(f"Wrote {CSV.name}: {len(assigned)} MAGs")
print(f"  case 1 GTDB genus       (_m): {cc['1_gtdb_genus']}")
print(f"  case 2 inherited ASV id      : {cc['2_asv_inherit']}")
print(f"  case 3 lowest GTDB label (_m): {cc['3_lowest_gtdb']}")
print(f"  duplicate iterativeIDs       : {len(dups)}  {dups}")
