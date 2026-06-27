#!/usr/bin/env python3
"""Recompute the root-level abundance Spearman correlation / p-value matrices
(correlation_matrix_root.csv, pvalue_matrix_root.csv) using the CURRENT ASV
iterativeID roots.  These two orphan files carried the old (numeric-root)
aggregation and had no generator left in the repo; this restores them on the
current nomenclature.  Roots = ASV iterativeID prefix (split('.')[0]); organisms
kept if their root-summed relative abundance reaches >= 1% in any sample (the
``_min_max_pct = 0.01`` filter from data_processing.py's root cell).

Run:  ~/Documents/py_venv/bin/python scripts/recompute_root_correlation_matrices.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
ab = pd.read_csv(ROOT / "abundances.csv").set_index("sample")
ab.columns = [c.split(".")[0] for c in ab.columns]            # current iterativeID roots
root_ab = ab.T.groupby(level=0).sum().T                       # samples x root
root_ab = root_ab.loc[:, root_ab.max(axis=0) >= 0.01]         # >=1% max-abundance filter
roots = sorted(root_ab.columns)
root_ab = root_ab[roots]

n = len(roots)
rho = np.ones((n, n)); pval = np.ones((n, n))                 # diagonal: rho=1, p=1 (self)
for i in range(n):
    for j in range(i + 1, n):
        r, p = spearmanr(root_ab.iloc[:, i], root_ab.iloc[:, j])
        if np.isnan(r):
            r, p = 0.0, 1.0
        rho[i, j] = rho[j, i] = r
        pval[i, j] = pval[j, i] = p

pd.DataFrame(rho, index=roots, columns=roots).to_csv(ROOT / "correlation_matrix_root.csv")
pd.DataFrame(pval, index=roots, columns=roots).to_csv(ROOT / "pvalue_matrix_root.csv")
print(f"wrote correlation_matrix_root.csv / pvalue_matrix_root.csv  ({n} roots, >=1% max)")
import re
print("numeric (old-style) roots remaining:", sum(1 for r in roots if re.search(r"\d", r)))
print("sample roots:", roots[:8])
