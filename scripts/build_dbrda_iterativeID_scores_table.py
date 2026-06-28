#!/usr/bin/env python3
"""Build a per-iterativeID table of the PERFORMANCE dbRDA ordination scores plus
the vector-based associations with each reactor-performance variable.

Sources (regenerate first with dbRDA_analysis.py + dbRDA_associations.py):
  - dbRDA/dbRDA_performance_Xs_5%_species.csv  -> ordination scores (plot_CAP1/2)
  - dbRDA/dbRDA_evidence_genus_vector_5%.csv   -> per (organism, vector) projection,
       cosine, |CAP|, rank, best-driver flag (panel == 'performance')

Output (dbRDA/):
  dbRDA_performance_iterativeID_scores.csv / .md
Each row = one iterativeID root (the dbRDA aggregates ASVs to their iterativeID
root); columns = ordination scores, peak phase, community module, the single
best-aligned performance driver, and the signed projection onto every
performance vector (association strength; sign = direction of alignment).
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent / "dbRDA"
sp = pd.read_csv(ROOT / "dbRDA_performance_Xs_5%_species.csv")
ev = pd.read_csv(ROOT / "dbRDA_evidence_genus_vector_5%.csv")
ev = ev[ev["panel"].str.contains("perf", case=False)].copy()

# vector order: keep the file's order, but surface the two new P-uptake dims
VEC_ORDER = ["specific denitrification rate", "denitrifying P-uptake rate",
             "anoxic:aerobic P ratio", "P removal", "N removal", "peak N2O"]
vecs = [v for v in VEC_ORDER if v in set(ev["vector"])] + \
       [v for v in sorted(ev["vector"].unique()) if v not in VEC_ORDER]

proj = ev.pivot_table(index="organism", columns="vector", values="projection", aggfunc="first")
proj = proj.reindex(columns=vecs)
proj.columns = [f"proj:{c}" for c in proj.columns]

meta = ev.groupby("organism").agg(
    module=("module", "first"),
    peak_phase=("peak_phase", "first"),
).reset_index()
best = (ev.sort_values("abs_CAP", ascending=False)
          .drop_duplicates("organism")[["organism", "vector"]]
          .rename(columns={"vector": "best_performance_driver"}))

tab = (sp.rename(columns={"genus": "iterativeID", "plot_CAP1": "CAP1", "plot_CAP2": "CAP2"})
         [["iterativeID", "functional_category", "CAP1", "CAP2"]]
         .merge(meta, left_on="iterativeID", right_on="organism", how="left").drop(columns="organism")
         .merge(best, left_on="iterativeID", right_on="organism", how="left").drop(columns="organism")
         .merge(proj, left_on="iterativeID", right_index=True, how="left"))

tab["_mag"] = tab[[c for c in tab.columns if c.startswith("proj:")]].abs().max(axis=1)
tab = tab.sort_values(["functional_category", "_mag"], ascending=[True, False]).drop(columns="_mag")
for c in tab.columns:
    if c.startswith("proj:") or c in ("CAP1", "CAP2"):
        tab[c] = tab[c].round(3)

out_csv = ROOT / "dbRDA_performance_iterativeID_scores.csv"
tab.to_csv(out_csv, index=False)
# compact markdown (top 40 by association magnitude)
md_cols = ["iterativeID", "functional_category", "CAP1", "CAP2", "peak_phase",
           "best_performance_driver"] + [c for c in tab.columns if c.startswith("proj:")]
top = tab.reindex(tab[[c for c in tab.columns if c.startswith('proj:')]].abs().max(axis=1)
                  .sort_values(ascending=False).index).head(40)[md_cols]
def _md(df):
    hdr = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
            for r in df.itertuples(index=False)]
    return "\n".join([hdr, sep] + rows) + "\n"
(ROOT / "dbRDA_performance_iterativeID_scores.md").write_text(_md(top))
print(f"wrote {out_csv.name} ({len(tab)} iterativeIDs x {len([c for c in tab.columns if c.startswith('proj:')])} perf vectors) + .md (top 40)")
print("performance vectors:", vecs)
print(tab.head(6).to_string(index=False))
