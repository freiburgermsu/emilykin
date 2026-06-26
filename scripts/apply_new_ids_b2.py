#!/usr/bin/env python3
"""Stage B2: relabel the remaining full-ASV-id-keyed / MAG-keyed artifacts that
Stage B did not cover (correlations/, network membership, allbins nosZ master).
All are bijective relabels (values unchanged).  ROOT-aggregated *_genera tables
are intentionally NOT touched here (they need re-aggregation, handled separately).
"""
import json, csv, glob, os, re
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
H = json.loads((Path(__file__).resolve().parent / "_new_ids.json").read_text())
a = H["asv_old2new"]; mag_new = H["mag_new"]
def rl(x): return a.get(x, x)
rep = []

# --- correlation matrices: header columns are full ASV ids -----------------
for name in ["correlation_n_matrix.csv","correlation_pvalue_matrix.csv",
             "correlation_qvalue_matrix.csv","correlation_rho_matrix.csv"]:
    p = ROOT/"correlations"/name
    if not p.exists(): continue
    rows=list(csv.reader(open(p)))
    n=sum(1 for c in rows[0] if c in a)
    rows[0]=[rl(c) for c in rows[0]]
    csv.writer(open(p,"w",newline="")).writerows(rows)
    rep.append(f"{name}: {n} header cols relabeled")

# --- ASV-column tables (NOT the _genera root-aggregated ones) --------------
for name in ["correlation_confirmed_dual.csv","correlation_dual_qvalue_table.csv",
             "correlation_dual_qvalue_table_extended.csv","correlation_long_table.csv",
             "correlation_robust_final.csv","correlation_robust_n20_q05.csv",
             "correlation_with_time_partial.csv","correlation_with_time_partial_extended.csv"]:
    p = ROOT/"correlations"/name
    if not p.exists(): continue
    rows=list(csv.reader(open(p))); hdr=rows[0]
    if "ASV" not in hdr: rep.append(f"{name}: no ASV col, skipped"); continue
    j=hdr.index("ASV"); n=0
    for r in rows[1:]:
        if len(r)>j and r[j] in a: r[j]=a[r[j]]; n+=1
    csv.writer(open(p,"w",newline="")).writerows(rows)
    rep.append(f"{name}: {n} ASV-column cells relabeled")

# --- per_parameter/*.json: top-level keys are full ASV ids -----------------
npp=0
for p in glob.glob(str(ROOT/"correlations"/"per_parameter"/"*.json")):
    d=json.load(open(p))
    if isinstance(d,dict):
        d={rl(k):v for k,v in d.items()}
        json.dump(d,open(p,"w")); npp+=1
rep.append(f"per_parameter/: {npp} json files relabeled")

# --- network/network_module_membership_*.json: member lists ----------------
nm=0
for p in glob.glob(str(ROOT/"network"/"network_module_membership_*.json")):
    d=json.load(open(p)); ch=False
    for mod,info in d.items():
        if isinstance(info,dict) and isinstance(info.get("members"),list):
            info["members"]=[rl(x) for x in info["members"]]; ch=True
    if ch: json.dump(d,open(p,"w"),indent=2); nm+=1
rep.append(f"network membership: {nm} json files relabeled")

# --- clade_classify/out/allbins_nosz_master.tsv: iterativeID col by MAG -----
p=ROOT/"clade_classify"/"out"/"allbins_nosz_master.tsv"
if p.exists():
    rows=list(csv.reader(open(p),delimiter="\t")); hdr=rows[0]
    mc=hdr.index("MAG"); ic=hdr.index("iterativeID"); n=0
    for r in rows[1:]:
        if len(r)>ic and r[mc] in mag_new: r[ic]=mag_new[r[mc]]; n+=1
    csv.writer(open(p,"w",newline=""),delimiter="\t").writerows(rows)
    rep.append(f"allbins_nosz_master.tsv: {n} iterativeID cells relabeled (MAG-keyed)")

print("\n".join(rep))
