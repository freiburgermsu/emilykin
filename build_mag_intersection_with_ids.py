#!/usr/bin/env python3
"""
Build a permutation of mag_abundance_by_day_intersection.csv that adds:
  - iterativeID (first column): either an existing ASV iterativeID whose
    taxonomy matches the MAG's GTDB taxonomy, or a new MAG-derived
    iterativeID (same prefix convention, _m suffix).
  - GTDB taxonomy columns (Kingdom–Species) — already present in the
    intersection CSV but ensured here from the canonical GTDB source.

iterativeID assignment rules (per MAG):
  1. Collect all primary-alignment ASVs that mapped to this MAG's 16S gene(s).
  2. For each ASV, retrieve its iterativeID and taxonomy.
  3. Check concordance with the MAG's GTDB taxonomy, evaluating from
     most-specific to least-specific rank (Species → Genus → Family →
     Order → Class → Phylum → Kingdom).
  4. Among concordant ASVs, pick the one whose iterativeID prefix is most
     specific (deepest matching rank). Ties broken by choosing the lowest
     numeric suffix so the assignment is deterministic.
  5. If no concordant ASV exists: mint a new ID as <prefix>.<N>_m where
     prefix = most-specific non-empty GTDB rank value, N = next unused
     integer for that prefix among all _m IDs generated this run.

Output: mag_abundance_by_day_intersection_with_ids.csv
"""
from collections import defaultdict
from pathlib import Path
import json
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
EMILYKIN   = Path("/scratch1/afreiburger/emilykin/EmilyKin")
ASV_MAG    = Path("/scratch1/afreiburger/emilykin/asv_mag_mapping")

INTERSECT  = EMILYKIN / "mag_abundance_by_day_intersection.csv"
IDS_JSON   = EMILYKIN / "iterativeIDs.json"
TAX_JSON   = EMILYKIN / "iterativeID_taxonomy.json"
SAM        = ASV_MAG  / "asv_vs_mag16s.sam"
OUT        = EMILYKIN / "mag_abundance_by_day_intersection_with_ids.csv"

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# ── 1. Load abundance + taxonomy table ────────────────────────────────────────
df = pd.read_csv(INTERSECT, index_col=0)
print(f"Intersection table: {len(df)} MAGs")

# Ensure all rank columns exist (fill missing with "")
for r in RANKS:
    if r not in df.columns:
        df[r] = ""
    df[r] = df[r].fillna("").astype(str)

# ── 2. Build ASV hash → iterativeID reverse map ───────────────────────────────
raw_ids  = json.load(open(IDS_JSON))          # {iterativeID: hash}
hash2id  = {v: k for k, v in raw_ids.items()} # {hash: iterativeID}

asv_tax  = json.load(open(TAX_JSON))          # {iterativeID: {Kingdom, ...}}
print(f"ASV iterativeIDs: {len(hash2id)}")

# ── 3. Parse SAM → {mag_name: set of primary-alignment ASV hashes} ────────────
# Primary = flag & 2304 == 0  (not secondary 256, not supplementary 2048)
mag_to_asvs: dict[str, set] = defaultdict(set)
with open(SAM) as fh:
    for line in fh:
        if line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        qname, flag, rname = parts[0], int(parts[1]), parts[2]
        if flag & 2304:          # secondary or supplementary → skip
            continue
        if flag & 4:             # unmapped → skip
            continue
        if rname == "*":
            continue
        # gene name format: CAN_1_bin.148__GDJGFH_00714
        mag = rname.split("__")[0] if "__" in rname else rname
        mag_to_asvs[mag].add(qname)

print(f"MAGs with mapped ASVs: {len(mag_to_asvs)}")

# ── 4. Helpers ────────────────────────────────────────────────────────────────

def _val(tax_dict: dict, rank: str) -> str:
    """Return cleaned taxonomy value; treat 'nan' as empty."""
    v = str(tax_dict.get(rank, "")).strip()
    return "" if v.lower() == "nan" else v


def most_specific_rank(tax_dict: dict) -> str:
    """Return the value of the deepest non-empty rank."""
    for rank in reversed(RANKS):
        v = _val(tax_dict, rank)
        if v:
            return v
    return "Bacteria"  # fallback


def most_specific_depth(tax_dict: dict) -> int:
    """Return the index of the deepest non-empty rank (-1 if none)."""
    depth = -1
    for i, rank in enumerate(RANKS):
        if _val(tax_dict, rank):
            depth = i
    return depth

def concordance_depth(asv_tax_dict: dict, mag_tax_dict: dict) -> int:
    """
    Return the 0-based match depth (index into RANKS) if:
      (a) There are no contradictions at any rank (no rank where both ASV
          and MAG have different non-empty values).
      (b) The ASV's own deepest non-empty rank agrees with the MAG.
      (c) The ASV's deepest classification reaches at least as deep as
          the MAG's deepest classification (prevents a Kingdom-only ASV
          from matching a Genus-level MAG).
    Returns -1 if any condition fails.
    """
    asv_depth = most_specific_depth(asv_tax_dict)
    mag_depth = most_specific_depth(mag_tax_dict)

    if asv_depth == -1:
        return -1  # ASV has no taxonomy

    # ASV must be at least as specific as the MAG's deepest resolved rank
    if asv_depth < mag_depth:
        return -1

    # The MAG must agree at the ASV's own most-specific level
    asv_val = _val(asv_tax_dict, RANKS[asv_depth])
    mag_val = _val(mag_tax_dict, RANKS[asv_depth])
    if mag_val != asv_val:
        return -1

    # No contradictions at any broader rank
    for i in range(asv_depth + 1):
        av = _val(asv_tax_dict, RANKS[i])
        mv = _val(mag_tax_dict, RANKS[i])
        if av and mv and av != mv:
            return -1

    return asv_depth

def numeric_suffix(iterative_id: str) -> int:
    """Extract the trailing integer from an iterativeID like 'Genus.42'."""
    try:
        return int(iterative_id.rsplit(".", 1)[-1])
    except ValueError:
        return 0

# ── 5. Assign iterativeID per MAG ────────────────────────────────────────────
m_counters: dict[str, int] = defaultdict(int)  # prefix → next _m counter
assigned: dict[str, str] = {}

for mag in df.index:
    mag_tax = {r: _val({r: df.at[mag, r]}, r) for r in RANKS}

    # Collect concordant ASV iterativeIDs
    best_depth  = -1
    best_ids    = []   # (numeric_suffix, iterativeID) tuples at best_depth

    for asv_hash in mag_to_asvs.get(mag, set()):
        iid = hash2id.get(asv_hash)
        if iid is None:
            continue
        atax = asv_tax.get(iid, {})
        depth = concordance_depth(atax, mag_tax)
        if depth > best_depth:
            best_depth = depth
            best_ids   = [(numeric_suffix(iid), iid)]
        elif depth == best_depth and depth >= 0:
            best_ids.append((numeric_suffix(iid), iid))

    if best_depth >= 0:
        # Use the concordant ASV iterativeID with the lowest numeric suffix
        best_ids.sort()
        assigned[mag] = best_ids[0][1]
    else:
        # Mint a new _m ID
        prefix = most_specific_rank(mag_tax)
        m_counters[prefix] += 1
        assigned[mag] = f"{prefix}.{m_counters[prefix]}_m"

# ── 6. Assemble output ────────────────────────────────────────────────────────
# Column order: iterativeID, day abundance columns, then GTDB ranks
day_cols  = [c for c in df.columns if c.startswith("day_")]
tax_cols  = RANKS

out = df[day_cols + tax_cols].copy()
out.insert(0, "iterativeID", [assigned[m] for m in out.index])

out.to_csv(OUT)
print(f"\nWritten: {OUT}  ({len(out)} rows × {len(out.columns)} columns)")

# ── 7. Summary stats ──────────────────────────────────────────────────────────
asv_reused = sum(1 for v in assigned.values() if not v.endswith("_m"))
mag_minted = sum(1 for v in assigned.values() if v.endswith("_m"))
print(f"  ASV iterativeID reused : {asv_reused}")
print(f"  New _m IDs minted      : {mag_minted}")

# Show depth distribution for reused IDs
depth_counts = defaultdict(int)
for mag, iid in assigned.items():
    if iid.endswith("_m"):
        continue
    mag_tax = {r: df.at[mag, r] for r in RANKS}
    atax = asv_tax.get(iid, {})
    d = concordance_depth(atax, mag_tax)
    depth_counts[RANKS[d] if d >= 0 else "none"] += 1

print(f"\n  Match depth distribution (reused ASV IDs):")
for rank in RANKS:
    if rank in depth_counts:
        print(f"    {rank:<12}: {depth_counts[rank]}")

print(f"\nSample assignments:")
for mag in list(df.index)[:10]:
    print(f"  {mag:<30} → {assigned[mag]}")
