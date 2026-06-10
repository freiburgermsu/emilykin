#!/usr/bin/env python3
"""
Build MAG quality metrics table from CheckM2 reports + GTDB taxonomy.

Rows = 276 dereplicated MAGs
Columns = Completeness, Contamination, Total_Contigs, Genome_Size_bp,
          Contig_N50, Max_Contig_Length, GC_Content, Coding_Density,
          Total_Coding_Sequences, Genome_Type, Domain, Phylum, Class,
          Order, Family, Genus, Species
"""
import os
import sys
from pathlib import Path
import pandas as pd

DREP_DIR   = Path("/scratch1/afreiburger/emilykin/processed/mag/drep/dereplicated_genomes")
CHECKM2    = Path("/scratch1/afreiburger/emilykin/processed/mag/checkm2")
GTDBTK_TSV = Path("/scratch1/afreiburger/emilykin/processed/mag/gtdbtk/classify/gtdbtk.bac120.summary.tsv")
OUT_TSV    = Path("/scratch1/afreiburger/emilykin/EmilyKin/mag_quality_summary.tsv")

# ── 1. Get the 276 dereplicated MAG names ──────────────────────────────────
drep_mags = sorted(p.stem for p in DREP_DIR.glob("*.fa"))
print(f"Dereplicated MAGs: {len(drep_mags)}")

# ── 2. Load CheckM2 reports; prefix names to match drep ───────────────────
# Per-sample: bin.X  → CAN_N_bin.X
# Coassembly:  bin.X → coasm_bin.X
sources = {
    "CAN_1": CHECKM2 / "per_sample_hybrid_CAN_1" / "quality_report.tsv",
    "CAN_2": CHECKM2 / "per_sample_hybrid_CAN_2" / "quality_report.tsv",
    "CAN_3": CHECKM2 / "per_sample_hybrid_CAN_3" / "quality_report.tsv",
    "CAN_4": CHECKM2 / "per_sample_hybrid_CAN_4" / "quality_report.tsv",
    "CAN_5": CHECKM2 / "per_sample_hybrid_CAN_5" / "quality_report.tsv",
    "coasm": CHECKM2 / "coassembly_hybrid"        / "quality_report.tsv",
}

frames = []
for prefix, path in sources.items():
    df = pd.read_csv(path, sep="\t")
    if prefix == "coasm":
        df["MAG"] = "coasm_" + df["Name"]
    else:
        df["MAG"] = prefix + "_" + df["Name"]
    frames.append(df)

checkm2 = pd.concat(frames, ignore_index=True)
checkm2 = checkm2.set_index("MAG")
print(f"CheckM2 total rows: {len(checkm2)}")

# Keep only dereplicated MAGs
missing_qc = [m for m in drep_mags if m not in checkm2.index]
if missing_qc:
    print(f"WARNING: {len(missing_qc)} derep MAGs not in CheckM2 reports: {missing_qc[:5]}...")
checkm2 = checkm2.loc[[m for m in drep_mags if m in checkm2.index]]
print(f"After filtering to drep set: {len(checkm2)}")

# ── 3. Load GTDB taxonomy ──────────────────────────────────────────────────
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PREFIX_MAP = dict(zip(["d__", "p__", "c__", "o__", "f__", "g__", "s__"], RANKS))

gtdb_raw = pd.read_csv(GTDBTK_TSV, sep="\t", usecols=["user_genome", "classification"])
gtdb_raw = gtdb_raw.set_index("user_genome")

def parse_gtdb(classif_str):
    parts = {v: "" for v in RANKS}
    for field in str(classif_str).split(";"):
        field = field.strip()
        for pfx, rank in PREFIX_MAP.items():
            if field.startswith(pfx):
                val = field[len(pfx):]
                parts[rank] = val if val else ""
    return pd.Series(parts)

gtdb_tax = gtdb_raw["classification"].apply(parse_gtdb)
print(f"GTDB rows: {len(gtdb_tax)}")

# ── 4. Select and rename CheckM2 columns ──────────────────────────────────
col_map = {
    "Completeness":             "Completeness_pct",
    "Contamination":            "Contamination_pct",
    "Completeness_Model_Used":  "Genome_Type",
    "Total_Contigs":            "Num_Contigs",
    "Genome_Size":              "Genome_Size_bp",
    "Contig_N50":               "N50_bp",
    "Max_Contig_Length":        "Max_Contig_Length_bp",
    "GC_Content":               "GC_Content",
    "Coding_Density":           "Coding_Density",
    "Total_Coding_Sequences":   "Total_CDS",
}
qc = checkm2[[c for c in col_map if c in checkm2.columns]].rename(columns=col_map)

# Simplify Genome_Type values
if "Genome_Type" in qc.columns:
    qc["Genome_Type"] = qc["Genome_Type"].str.replace(
        r"Neural Network \(Specific Model\)", "Bacteria", regex=True
    ).str.replace(
        r"Neural Network \(General Model\)", "Archaea", regex=True
    )

# ── 5. Join ────────────────────────────────────────────────────────────────
result = qc.join(gtdb_tax, how="left")
result.index.name = "MAG"

# Round float columns for readability
for col in ["Completeness_pct", "Contamination_pct", "GC_Content", "Coding_Density"]:
    if col in result.columns:
        result[col] = result[col].round(2)

# Sort by completeness desc
result = result.sort_values("Completeness_pct", ascending=False)

# ── 6. Write ───────────────────────────────────────────────────────────────
result.to_csv(OUT_TSV, sep="\t")
print(f"\nWritten: {OUT_TSV}  ({len(result)} MAGs × {len(result.columns)} columns)")
print(f"\nColumn list: {list(result.columns)}")
print(f"\nFirst 5 rows:")
print(result.head(5).to_string())
