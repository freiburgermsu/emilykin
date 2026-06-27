#!/usr/bin/env python3
"""
regen_iterative_ids.py  —  Stage A of the iterativeID nomenclature change.

NEW RULES
=========
ASV iterativeID root  (per ASV, from SILVA/MiDAS taxonomy in taxonomy.json):
    The root is the DEEPEST taxonomic level among Kingdom..Genus (Species is
    excluded, preserving the genus-level convention) whose value is
        * non-empty, AND
        * contains NO digit 0-9  (numeric placeholders such as midas_g_171,
          SM1A02, Pir4_lineage, Christensenellaceae_R-7_group are skipped and
          the next-shallower clean label is used).
    Fallback = the Kingdom value ('Bacteria'/'Archaea') or 'Bacteria'.
    iterativeID = f"{root}.{N}"  (N = per-root running count, ASV table order).

MAG iterativeID  (per MAG, GTDB taxonomy in mag_abundance_summary.json):
    case 1  GTDB Genus defined & non-numeric          -> f"{Genus}.{N}_m"
    case 2  else, a mapped (non-weak) ASV agrees with the MAG's GTDB taxonomy
            at the DEEPEST defined non-numeric GTDB rank (Phylum..Genus; Domain
            is never a valid match rank); pick the agreeing ASV with the highest
            mean relative abundance -> the ASV's own new iterativeID (NO suffix)
    case 3  else (genus not usable AND no agreeing ASV)
            lowest defined non-numeric GTDB label      -> f"{label}.{N}_m"
    '_m' is appended AFTER the per-root running count; inherited IDs (case 2)
    carry no suffix (the MAG and its representative ASV share a label).

Outputs (repo root):
    iterativeIDs.json, iterativeID_levels.json, iterativeID_phylums.json,
    iterativeID_taxonomy.json, taxonIDs.json, iterativeID_color_map.json
    asv_iterativeID_old_to_new.csv / .json
    mag_iterativeID_old_to_new.csv / .json
and writes scripts/_new_ids.json (machine-readable handoff for Stage B):
    {"asv_old2new":{...}, "asv_hash2new":{...}, "mag_new":{mag:new_id},
     "mag_case":{mag:case}, "asv_new2hash":{...}}

Run:  ~/Documents/py_venv/bin/python scripts/regen_iterative_ids.py
"""
from collections import defaultdict, Counter
from pathlib import Path
import csv, json, re

ROOT = Path(__file__).resolve().parent.parent
HANDOFF = Path(__file__).resolve().parent / "_new_ids.json"

HASDIGIT = re.compile(r"\d")
def numeric(s) -> bool:
    return bool(HASDIGIT.search(str(s)))
def clean(v) -> str:
    s = str(v).strip()
    return "" if s.lower() in ("", "nan", "none") else s
def norm(v) -> str:
    # Compensate for naming-convention differences (Ca_ prefixes, GTDB _A/_F/_m
    # suffixes, sensu_stricto, ...) by comparing the longest token of split('_').
    v = clean(v)
    return max(v.split("_"), key=len).lower() if v else ""

ASV_LEVELS  = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]  # no Species
GTDB_LEVELS = ["Domain",  "Phylum", "Class", "Order", "Family", "Genus"]  # no Species
# GTDB 'Domain' corresponds to the ASV (SILVA) 'Kingdom' rank when matching.
GTDB2ASV_RANK = {"Domain": "Kingdom", "Phylum": "Phylum", "Class": "Class",
                 "Order": "Order", "Family": "Family", "Genus": "Genus"}

# ── load sources ──────────────────────────────────────────────────────────────
taxonomy = json.loads((ROOT / "taxonomy.json").read_text())            # hash -> {rank:val}
old_ids  = json.loads((ROOT / "iterativeIDs.json").read_text())        # old_iid -> hash
relab    = json.loads((ROOT / "relative_abundance.json").read_text())  # hash -> {sample:relab}
a2m      = json.loads((ROOT / "asv_to_mag_mapping.json").read_text())["asv_to_mag"]
mags     = json.loads((ROOT / "mag_abundance_summary.json").read_text())  # mag -> {taxonomy,...}

mean_ab  = {h: (sum(d.values()) / len(d) if d else 0.0) for h, d in relab.items()}
old_hash2id = {h: i for i, h in old_ids.items()}

# canonical MAG order = mag_abundance_by_day_intersection_with_ids.csv row order
with open(ROOT / "mag_abundance_by_day_intersection_with_ids.csv") as fh:
    mag_order = [r[0] for r in csv.reader(fh)][1:]
mag_order = [m for m in mag_order if m]

# old canonical MAG ids (mag_abundance_summary.tsv 'mag_iterativeID' col) for the mapping
old_mag_id = {}
old_mag_intersection_id = {}
with open(ROOT / "mag_abundance_summary.tsv") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        old_mag_id[r["MAG"]] = r.get("mag_iterativeID", "")
with open(ROOT / "mag_abundance_by_day_intersection_with_ids.csv") as fh:
    rd = csv.DictReader(fh)
    for r in rd:
        old_mag_intersection_id[r["MAG"]] = r.get("iterativeID", "")

# ── ASV iterativeIDs (new) ────────────────────────────────────────────────────
def asv_root_level(tax):
    chosen = None
    for L in ASV_LEVELS:
        v = clean(tax.get(L, ""))
        if v and not numeric(v):
            chosen = (v, L)
    if chosen is None:
        v = clean(tax.get("Kingdom", "")) or "Bacteria"
        chosen = (v, "Kingdom")
    return chosen  # (root_value, level_name)

new_iterativeIDs   = {}   # new_iid -> hash
new_levels         = {}   # new_iid -> level
new_phylums        = {}   # new_iid -> Phylum value
new_taxonomy       = {}   # new_iid -> {rank:val}
new_taxonIDs       = defaultdict(list)   # root -> [hash,...]
hash2new           = {}   # hash -> new_iid
root_counter       = defaultdict(int)

for h, tax in taxonomy.items():
    root, lvl = asv_root_level(tax)
    root_counter[root] += 1
    iid = f"{root}.{root_counter[root]}"
    new_iterativeIDs[iid] = h
    new_levels[iid]       = lvl
    new_phylums[iid]      = tax.get("Phylum", "")
    new_taxonomy[iid]     = tax
    new_taxonIDs[root].append(h)
    hash2new[h]           = iid

# ASV old -> new (join by hash)
asv_old2new = {}
for old_iid, h in old_ids.items():
    asv_old2new[old_iid] = hash2new.get(h, "")

# ── MAG iterativeIDs (new) ────────────────────────────────────────────────────
NONWEAK = {"species", "genus", "family"}
mag_asvs = defaultdict(list)   # mag -> [hash,...]  (non-weak mapped ASVs)
for h, rec in a2m.items():
    if rec.get("confidence_tier") in NONWEAK and rec.get("best_mag"):
        mag_asvs[rec["best_mag"]].append(h)

def gtdb(mag, rank):
    return clean(mags.get(mag, {}).get("taxonomy", {}).get(rank, ""))

def deepest_nonnum_gtdb(mag):
    chosen = ("", "")
    for L in GTDB_LEVELS:
        v = gtdb(mag, L)
        if v and not numeric(v):
            chosen = (L, v)
    return chosen   # (rank, value)

mag_new   = {}
mag_case  = {}
m_counter = defaultdict(int)

for mag in mag_order:
    # Priority: a mapped ASV that AGREES with the MAG's GTDB taxonomy (at the
    # deepest defined non-numeric GTDB rank) wins over the GTDB-derived label.
    r_rank, r_val = deepest_nonnum_gtdb(mag)
    agree = []
    if r_rank and r_rank != "Domain" and r_val:     # no trivial Domain match
        asv_rank = GTDB2ASV_RANK[r_rank]
        for h in mag_asvs.get(mag, []):
            av = clean(taxonomy.get(h, {}).get(asv_rank, ""))
            if av and norm(av) == norm(r_val):
                agree.append(h)
    if agree:                                        # >=1 mapped ASV agrees -> inherit
        best = max(agree, key=lambda h: (mean_ab.get(h, 0.0), hash2new.get(h, "")))
        mag_new[mag]  = hash2new[best]               # most-abundant AGREEING ASV id, no suffix
        mag_case[mag] = "asv_inherit"
        continue
    # No agreeing ASV (all mapped ASVs disagree, or none mapped) -> GTDB-defined label
    g = gtdb(mag, "Genus")
    if g and not numeric(g):                         # GTDB genus defined & non-numeric
        m_counter[g] += 1
        mag_new[mag]  = f"{g}.{m_counter[g]}_m"
        mag_case[mag] = "gtdb_genus"
    else:                                            # lowest defined non-numeric GTDB label
        lbl = r_val if r_val else "Bacteria"
        m_counter[lbl] += 1
        mag_new[mag]  = f"{lbl}.{m_counter[lbl]}_m"
        mag_case[mag] = "gtdb_lowest"

# ── colour map (replicates data_processing.py cell @169-210 with new ids) ──────
import matplotlib.pyplot as plt
def rgba_to_hex(rgba):
    return "#{:02x}{:02x}{:02x}".format(*(int(c * 255) for c in rgba[:3]))
phyl = new_phylums
archaea_phyla = sorted({v for v in phyl.values()
                        if v and ("archaeo" in v.lower()
                        or any(a in v.lower() for a in
                               ["candidatus thermoplasmatota", "halobacterota", "methanobacteriota"]))})
bacteria_phyla = sorted({v for v in phyl.values() if v and v not in archaea_phyla})
taxa_color_map = {}
for i, p in enumerate(archaea_phyla):
    taxa_color_map[p] = rgba_to_hex(plt.cm.turbo(i / max(len(archaea_phyla), 1) * 0.15))
for i, p in enumerate(bacteria_phyla):
    taxa_color_map[p] = rgba_to_hex(plt.cm.turbo(0.2 + i / max(len(bacteria_phyla), 1) * 0.8))
iterativeID_color_map = {iid: taxa_color_map[p] for iid, p in phyl.items() if p}
# expand Proteobacteria into per-class Purples
def expand(phylum, cmap, lo=0.6, hi=0.95, base_t=0.85):
    classes = sorted({t["Class"] for t in new_taxonomy.values()
                      if t.get("Phylum") == phylum and t.get("Class")})
    n = max(len(classes), 1)
    class_colors = {c: rgba_to_hex(cmap(lo + (hi - lo) * i / max(n - 1, 1))) for i, c in enumerate(classes)}
    for iid, t in new_taxonomy.items():
        if t.get("Phylum") == phylum and t.get("Class") in class_colors:
            iterativeID_color_map[iid] = class_colors[t["Class"]]
    taxa_color_map[phylum] = rgba_to_hex(cmap(base_t))
    return class_colors, rgba_to_hex(cmap(base_t))
proteo_class_color, proteo_base = expand("Proteobacteria", plt.cm.Purples)

# ── write outputs ─────────────────────────────────────────────────────────────
def dump(obj, name):
    (ROOT / name).write_text(json.dumps(obj, indent=0 if name.endswith("taxonomy.json") else None))

(ROOT / "iterativeIDs.json").write_text(json.dumps(new_iterativeIDs))
(ROOT / "iterativeID_levels.json").write_text(json.dumps(new_levels))
(ROOT / "iterativeID_phylums.json").write_text(json.dumps(new_phylums))
(ROOT / "iterativeID_taxonomy.json").write_text(json.dumps(new_taxonomy))
(ROOT / "taxonIDs.json").write_text(json.dumps({k: v for k, v in new_taxonIDs.items()}))
(ROOT / "iterativeID_color_map.json").write_text(json.dumps(iterativeID_color_map))
(ROOT / "proteo_class_color.json").write_text(json.dumps(proteo_class_color))
(ROOT / "phylum_base_overrides.json").write_text(json.dumps({"Proteobacteria": proteo_base}))
# Genus_color_map.json / Species_color_map.json = the phylum colour table (legacy names)
(ROOT / "Genus_color_map.json").write_text(json.dumps(taxa_color_map))
(ROOT / "Species_color_map.json").write_text(json.dumps(taxa_color_map))

# mapping files
with open(ROOT / "asv_iterativeID_old_to_new.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["asv_hash", "old_iterativeID", "new_iterativeID",
                                    "old_root", "new_root", "changed"])
    for old_iid, h in old_ids.items():
        new = asv_old2new[old_iid]
        oro, nro = old_iid.rsplit(".", 1)[0], new.rsplit(".", 1)[0]
        w.writerow([h, old_iid, new, oro, nro, "yes" if oro != nro else "no"])
(ROOT / "asv_iterativeID_old_to_new.json").write_text(json.dumps(asv_old2new))

with open(ROOT / "mag_iterativeID_old_to_new.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["MAG", "old_mag_iterativeID", "old_intersection_id",
                                    "new_mag_iterativeID", "assignment_case"])
    for mag in mag_order:
        w.writerow([mag, old_mag_id.get(mag, ""), old_mag_intersection_id.get(mag, ""),
                    mag_new[mag], mag_case[mag]])
(ROOT / "mag_iterativeID_old_to_new.json").write_text(
    json.dumps({mag: mag_new[mag] for mag in mag_order}))

HANDOFF.write_text(json.dumps({
    "asv_old2new": asv_old2new,
    "asv_hash2new": hash2new,
    "asv_new2hash": new_iterativeIDs,
    "mag_new": mag_new,
    "mag_case": mag_case,
}))

# ── verification report ───────────────────────────────────────────────────────
print("=== ASV ===")
print("  ASVs:", len(new_iterativeIDs), "| unique new iids:", len(set(new_iterativeIDs)),
      "| unique hashes:", len(set(new_iterativeIDs.values())))
asv_changed = sum(1 for o, n in asv_old2new.items() if o.rsplit('.',1)[0] != n.rsplit('.',1)[0])
print("  old->new bijective:", len(asv_old2new) == len(set(asv_old2new.values())),
      "| roots changed:", asv_changed)
numroots_new = [r for r in new_taxonIDs if numeric(r)]
print("  numeric NEW roots (should be 0):", len(numroots_new), numroots_new[:5])
print("  top new roots:", [f"{r}:{len(new_taxonIDs[r])}" for r, _ in
                           Counter({k: len(v) for k, v in new_taxonIDs.items()}).most_common(8)])
print("=== MAG ===")
cc = Counter(mag_case.values())
print("  cases:", dict(cc), "| total:", len(mag_new))
dups = {v: c for v, c in Counter(mag_new.values()).items() if c > 1}
print("  duplicate MAG ids:", len(dups))
for v, c in list(dups.items())[:20]:
    holders = [m for m in mag_order if mag_new[m] == v]
    print(f"      {v} x{c}  {holders}")
print("  numeric-root _m ids (should be 0):",
      sum(1 for m in mag_order if mag_new[m].endswith("_m") and numeric(mag_new[m].rsplit('.',1)[0])))
print("  sample MAG assignments:")
for m in mag_order[:12]:
    print(f"      {m:18s} {old_mag_id.get(m,''):24s} -> {mag_new[m]:24s} [{mag_case[m]}]")
