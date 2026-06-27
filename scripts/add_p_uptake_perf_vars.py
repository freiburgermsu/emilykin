#!/usr/bin/env python3
"""Add two batch-assay P-uptake performance dimensions to
performance_data_with_sample_ids.csv for use as dbRDA performance constraint
vectors.  Values are the PER-PHASE figures the researcher curated on worksheet 2
('addition phase perfromance data') of 'Performance data_Andrew.xlsx':

  1. specific denitrifying P uptake rate   (max-range / first-hour anoxic uptake)
  2. anoxic:aerobic P uptake rate ratio     (anoxic vs aerobic P-uptake rate)

Each dbRDA sample inherits its phase's value (phases I-V via phase_day_ranges.json).
Phase III has no anoxic:aerobic ratio reported; those samples are left blank and
the dbRDA fills them by phase-median (overall-median fallback) imputation.

Run:  ~/Documents/py_venv/bin/python scripts/add_p_uptake_perf_vars.py
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import openpyxl

REPO = Path(__file__).resolve().parent.parent
CSV = REPO / "performance_data_with_sample_ids.csv"
XLSX = REPO / "Performance data_Andrew.xlsx"
V1 = "specific denitrifying P uptake rate"
V2 = "anoxic:aerobic P uptake rate ratio"

# --- per-phase values from worksheet 2 ---
wb = openpyxl.load_workbook(XLSX, data_only=True)
ws = wb["addition phase perfromance data"]
per_phase = {}
for ph, v1, v2, *_ in ws.iter_rows(min_row=2, values_only=True):
    if ph and str(ph).strip() in ("I", "II", "III", "IV", "V"):
        per_phase[str(ph).strip()] = (
            pd.to_numeric(v1, errors="coerce"),
            pd.to_numeric(v2, errors="coerce"),
        )
print("per-phase values (phase: v1, v2):", per_phase)

# --- map each performance row to its phase ---
phr = json.loads((REPO / "phase_day_ranges.json").read_text())
def day_to_phase(d):
    if pd.isna(d):
        return None
    for ph, span in phr.items():
        if span["start"] <= d <= span["end"]:
            return ph
    return None

df = pd.read_csv(CSV, low_memory=False)
df = df.drop(columns=[c for c in df.columns
                      if c == V1 or c == V2 or c.startswith("specific denitrifying P uptake")])
phase = pd.to_numeric(df["abundance_day"], errors="coerce").map(day_to_phase)
df[V1] = phase.map(lambda p: per_phase.get(p, (np.nan, np.nan))[0])
df[V2] = phase.map(lambda p: per_phase.get(p, (np.nan, np.nan))[1])
df.to_csv(CSV, index=False)

print(f"\nwrote {V1!r} and {V2!r} to {CSV.name}")
on_sid = df[df["sample_id"].notna()]
for c in (V1, V2):
    print(f"  {c}: non-null on dbRDA sample rows = {on_sid[c].notna().sum()}/{len(on_sid)} "
          f"| distinct phase values = {sorted(df[c].dropna().unique())}")
