#!/usr/bin/env python3
"""
assign_mag_iterative_ids.py

Reassign the ``iterativeID`` column of
``mag_abundance_by_day_intersection_with_ids.csv`` (one iterativeID per MAG).

Assignment rules (in priority order):

  Option 1  — inherit a mapped ASV's iterativeID   (suffix ``_ma``)
      a) Consider the ASVs that map to this MAG at species/genus/family
         confidence tier (asv_to_mag_mapping.json; the ``weak`` tier is
         excluded — this is exactly the set in the ``iterativeIDs`` column of
         mag_abundance_summary.tsv).
      b) Keep those whose ASV *Family* equals the MAG's GTDB *Family*
         (case-insensitive) — i.e. the ASV taxonomy matches the GTDB taxonomy
         to at least the Family level.
      c) Of the survivors, choose the ASV with the highest mean relative
         abundance across all samples (relative_abundance.json).
      => assigned ID = ``<ASV_iterativeID>_ma``

  Option 2  — mint a MAG-derived iterativeID        (suffix ``_m``)
      Used when the MAG has no mapped ASV, or none of the mapped ASVs agree at
      the Family level (covers novel GTDB lineages — UBA*, JAEDAM*, ... that
      have no SILVA/MiDAS ASV counterpart).
      prefix = deepest non-empty GTDB rank among
               Genus -> Family -> Order -> Class -> Phylum -> Domain
               (Species is skipped so the prefix stays a clean genus/family name),
      numbered per-prefix in MAG row order.
      => assigned ID = ``<prefix>.<N>_m``

GTDB taxonomy is read from the authoritative ``mag_abundance_summary.json``.
The Kingdom..Species columns *inside* the CSV are sparse and are NOT used for
matching; they are left untouched.  Only the ``iterativeID`` column is rewritten
(every other field is preserved verbatim), so the diff is limited to that column.

Run:  ~/Documents/py_venv/bin/python assign_mag_iterative_ids.py
"""
from collections import defaultdict, Counter
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent

CSV       = ROOT / "mag_abundance_by_day_intersection_with_ids.csv"
IDS_JSON  = ROOT / "iterativeIDs.json"            # {iterativeID: asv_hash}
TAX_JSON  = ROOT / "iterativeID_taxonomy.json"    # {iterativeID: {rank: value}}
ASV_MAG   = ROOT / "asv_to_mag_mapping.json"      # {asv_hash: {best_mag, confidence_tier, ...}}
MAG_SUM   = ROOT / "mag_abundance_summary.json"   # {mag: {taxonomy: {...}, abundance: {...}}}
RELAB     = ROOT / "relative_abundance.json"      # {asv_hash: {sample: rel_abundance}}

NONWEAK   = {"species", "genus", "family"}
# Ranks used to pick an Option-2 prefix, deepest first.  GTDB summary uses
# "Domain" (not "Kingdom"); Species is intentionally omitted.
PREFIX_RANKS = ["Genus", "Family", "Order", "Class", "Phylum", "Domain"]


def _clean(v) -> str:
    """Strip a taxonomy value; treat missing / 'nan' / 'none' as empty."""
    s = str(v).strip()
    return "" if s.lower() in ("", "nan", "none") else s


# ── load sources ──────────────────────────────────────────────────────────────
ids     = json.loads(IDS_JSON.read_text())          # iid -> hash
hash2id = {h: i for i, h in ids.items()}            # hash -> iid
itax    = json.loads(TAX_JSON.read_text())          # iid -> {rank: value}
a2m     = json.loads(ASV_MAG.read_text())["asv_to_mag"]
mags    = json.loads(MAG_SUM.read_text())           # mag -> {taxonomy, abundance}
relab   = json.loads(RELAB.read_text())             # hash -> {sample: relab}

# mean relative abundance per ASV hash (across all samples)
mean_ab = {h: (sum(d.values()) / len(d) if d else 0.0) for h, d in relab.items()}

# mapped ASVs per MAG, non-weak tiers only, as iterativeIDs
mag_to_iids: dict[str, list[str]] = defaultdict(list)
for h, rec in a2m.items():
    if rec.get("confidence_tier") not in NONWEAK:
        continue
    mag = rec.get("best_mag")
    iid = hash2id.get(h)
    if mag and iid:
        mag_to_iids[mag].append(iid)


def gtdb(mag: str, rank: str) -> str:
    return _clean(mags.get(mag, {}).get("taxonomy", {}).get(rank, ""))


def mag_prefix(mag: str) -> str:
    for r in PREFIX_RANKS:
        v = gtdb(mag, r)
        if v:
            return v
    return "Bacteria"


def assign(mag: str, m_counter: dict[str, int]) -> tuple[str, str]:
    """Return (iterativeID, kind) where kind is 'ma' or 'm'."""
    mag_fam = gtdb(mag, "Family")
    candidates = []  # (mean_abundance, iterativeID)
    if mag_fam:
        for iid in mag_to_iids.get(mag, []):
            asv_fam = _clean(itax.get(iid, {}).get("Family", ""))
            if asv_fam and asv_fam.lower() == mag_fam.lower():     # Family-level match
                candidates.append((mean_ab.get(ids.get(iid), 0.0), iid))
    if candidates:
        # most abundant; deterministic tie-break by iterativeID name
        candidates.sort(key=lambda t: (-t[0], t[1]))
        return f"{candidates[0][1]}_ma", "ma"
    prefix = mag_prefix(mag)
    m_counter[prefix] += 1
    return f"{prefix}.{m_counter[prefix]}_m", "m"


# ── compute assignments in CSV row order ──────────────────────────────────────
lines = CSV.read_text().splitlines(keepends=True)
header, body = lines[0], lines[1:]

assigned: dict[str, str] = {}
kinds: dict[str, str] = {}
m_counter: dict[str, int] = defaultdict(int)
for line in body:
    mag = line.split(",", 1)[0]
    if not mag:
        continue
    iid, kind = assign(mag, m_counter)
    assigned[mag] = iid
    kinds[mag] = kind

# ── rewrite ONLY the iterativeID column (field index 1), verbatim otherwise ───
out_lines = [header]
for line in body:
    nl = "\n" if line.endswith("\n") else ""
    core = line[: -len(nl)] if nl else line
    if not core:
        out_lines.append(line)
        continue
    mag, _old, rest = core.split(",", 2)   # MAG names never contain commas
    out_lines.append(f"{mag},{assigned[mag]},{rest}{nl}")

CSV.write_text("".join(out_lines))

# ── summary ───────────────────────────────────────────────────────────────────
n_ma = sum(1 for k in kinds.values() if k == "ma")
n_m  = sum(1 for k in kinds.values() if k == "m")
dups = {v: c for v, c in Counter(assigned.values()).items() if c > 1}
print(f"Wrote {CSV.name}: {len(assigned)} MAGs")
print(f"  Option 1  inherited ASV iterativeID (_ma): {n_ma}")
print(f"  Option 2  MAG-derived iterativeID    (_m): {n_m}")
print(f"  duplicate iterativeIDs across MAGs       : {len(dups)}")
for v, c in dups.items():
    print(f"      {v}  x{c}")
