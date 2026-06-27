"""Organism<->vector association spreadsheets from the db-RDA results.

For each inclusion threshold and each panel (operational drivers / performance),
quantify how every organism associates with every constraint vector:

  projection = |CAP| x cos   (association STRENGTH = how far the organism
                              projects onto the arrow; signed)
  cosine     = directional alignment (the earlier metric; misleading alone for
                              near-origin organisms)
  |CAP|      = distance from the origin (small => weak / non-specific)
  proj_rank  = rank among organisms for that vector (1 = strongest)
  best_driver= the single vector each organism aligns with best (argmax cos)
  module     = Louvain module (clustered in CAP space)
  peak_phase = phase of maximum mean relative abundance

Thresholds:
  base  -> top-10-genera-per-phase union (29 organisms)   [suffix '']
  5%    -> every genus whose max per-sample rel. abundance >= 5%  [suffix '_5%']

Outputs (graphs/dbRDA/):
  dbRDA_evidence_genus_vector{suffix}.csv                 (tidy: organisms x vectors)
  dbRDA_environmental_vector_associations{suffix}.csv
  dbRDA_performance_vector_associations{suffix}.csv

Base association tables are curated (hand-written interpretations from the 29-genus
fit). The _5% tables carry those interpretations + main phases, and recompute the
strongest associated genera over the broader pool, marking with '*' every genus
that lies BEYOND the original top-10 union.

Run (after dbRDA_analysis.py):
    ~/Documents/py_venv/bin/python dbRDA_associations.py
"""
from __future__ import annotations
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from pathlib import Path
import dbRDA_analysis as M

OUT = M.OUT_DIR
TOPN_5PCT = 12          # max genera listed per vector in the _5% tables
COS_MIN = 0.5           # min directional alignment for a genus to be "associated"

# pretty genus names matching the manuscript tables
DISP = {'Ca_Accumulibacter': 'Ca. Accumulibacter', 'Ca_Competibacter': 'Ca. Competibacter',
        'Ca_Sarcinithrix': 'Ca. Sarcinithrix', 'Ca_Caldilinea': 'Ca. Caldilinea',
        'Ca_Obscuribacter': 'Ca. Obscuribacter', 'Ca_Phosphoribacter': 'Ca. Phosphoribacter',
        'alphaI_cluster': 'alpha I-cluster'}
def disp(g: str) -> str: return DISP.get(g, g)

# ---- curated base tables (29-genus fit), N_Ax-1 dropped; N/P kept ------------
CURATED_ENV = [
 ("DO avg", "I", "M1",
  "Higher oxygen exposure selected the aeration-associated early community."),
 ("Cumulative aeration", "I", "M1",
  "Aeration history explains the early community state; near-collinear with DO avg and N_Ax-2."),
 ("N_Ax-2", "I-II", "M1",
  "Points with the aeration/early cluster rather than a distinct anoxic group - a Phase I->II aerobic-to-anoxic transition signal."),
 ("Acetate", "II", "M1",
  "Acetate selected GAO-like (Ca. Competibacter) and denitrification/carbon-cycling taxa."),
 ("C/N", "III/II", "M1",
  "C/N ratio structured denitrification/carbon-storage taxa."),
 ("Propionate", "V", "M3",
  "Phase V propionate-responsive denitrifying/PAO community."),
 ("N/P", "IV (-V)", "M2/M3",
  "Nutrient stoichiometry structured the Phase IV-V nutrient-removal community, centered on the PAO Ca. Accumulibacter. (N_Ax-1, near-collinear with N/P, was dropped.)"),
]
# curated genera per env vector (base 29-genus fit), for the base table
CURATED_ENV_GENERA = {
 "DO avg": "Thauera; Kapabacteriales; Sediminibacterium; Zoogloea; Pseudoxanthomonas",
 "Cumulative aeration": "Thauera; Zoogloea; Kapabacteriales; Sediminibacterium; Pseudoxanthomonas",
 "N_Ax-2": "Pseudoxanthomonas; Thauera; Zoogloea; Kapabacteriales; Sediminibacterium; Ca. Competibacter",
 "Acetate": "Ca. Competibacter; Azospira; Hydrogenophaga; Dokdonella; Denitratisoma",
 "C/N": "Denitratisoma; Dokdonella; Azospira; Ca. Competibacter; Hydrogenophaga",
 "Propionate": "Rhodoplanes; Ca. Sarcinithrix; alpha I-cluster; Ca. Accumulibacter; Rhizobiales; Gracilibacteria; Rhodospirillales; Holophagae",
 "N/P": "Ca. Accumulibacter; Ignavibacteria; Ca. Caldilinea; Flexifilaceae; Lysobacter; Saprospiraceae",
}
CURATED_PERF = [
 ("Specific denitrification rate", "V", "M3",
  "Later low-aeration operation selected a Phase-V community with stronger denitrification kinetics."),
 ("N removal (concentration)", "IV-V", "M3/M2",
  "Absolute N removal in the later low-aeration phases tracked the Ca. Accumulibacter / Phase-V denitrifying community."),
 ("P removal", "IV", "M2/M3",
  "P removal centered on the PAO Ca. Accumulibacter plus the Phase IV nutrient-removal cluster."),
 ("Peak N2O", "III-IV", "M2 (M1 edge)",
  "N2O accumulation tracked the Phase III-IV incomplete-denitrification/nutrient-removal cluster. Aeromonas points purest toward N2O but is near-origin (weak). Azospira and Ca. Competibacter are NOT N2O-associated."),
]
CURATED_PERF_GENERA = {
 "Specific denitrification rate": "Rhodoplanes; Ca. Sarcinithrix; Ca. Accumulibacter; alpha I-cluster; Rhizobiales; Gracilibacteria; Rhodospirillales; Anaerolineae; Holophagae",
 "N removal (concentration)": "Ca. Accumulibacter; Gracilibacteria; Rhodoplanes; Ca. Sarcinithrix; Flexifilaceae; Ca. Caldilinea; Lysobacter; Ca. Falkowbacteria",
 "P removal": "Ca. Accumulibacter; Ignavibacteria; Gracilibacteria; Lysobacter; Flexifilaceae; Ca. Caldilinea",
 "Peak N2O": "Ignavibacteria; Ca. Accumulibacter; Lysobacter; Ca. Falkowbacteria; Gracilibacteria; Denitratisoma",
}
# map display vector name -> raw biplot column name (performance panel)
PERF_VEC_RAW = {
 "Specific denitrification rate": "specific denitrification rates",
 "N removal (concentration)": "N removal (ppm) [N-ppn]",
 "P removal": "P removal [P%]",
 "Peak N2O": "peakN2O [mg/L]",
}
# env biplot labels match the curated names except for the lower-case aeration one
ENV_VEC_RAW = {"Cumulative aeration": "cumulative aeration"}


def peak_phase_map(long_file: Path) -> dict:
    tab = pd.read_csv(long_file, low_memory=False)
    tab = tab[tab['Phase'].isin(['I', 'II', 'III', 'IV', 'V'])]
    tab = tab[tab['Genus'].notna()].copy()
    tab['root'] = tab['seq'].map(M.SEQ2ROOT)   # key by iterativeID root, matching the ordination
    tab = tab[tab['root'].notna()]
    gp = (tab.groupby(['root', 'Phase'])['rel_ab'].mean().reset_index())
    order = ['I', 'II', 'III', 'IV', 'V']
    out = {}
    for g, sub in gp.groupby('root'):
        s = sub.set_index('Phase')['rel_ab'].reindex(order).fillna(0)
        out[g] = s.idxmax() if s.sum() > 0 else ''
    return out


def build_panels(Y, sample_ids, stages):
    stage_of = dict(zip(sample_ids, stages))
    Xp = M.build_perf_X_per_sample(sample_ids, stages).rename(columns=M.PERF_LABELS)
    kp = Xp.notna().all(axis=1)
    resP = M.dbrda(Y.loc[kp], Xp.loc[kp]); stP = [stage_of[s] for s in Y.loc[kp].index]
    Xi = M.build_influence_X_per_sample(sample_ids, stages)
    ke = Xi.notna().all(axis=1)
    resE = M.dbrda(Y.loc[ke], Xi.loc[ke]); stE = [stage_of[s] for s in Y.loc[ke].index]
    return ('environment', resE, stE), ('performance', resP, stP)


def panel_metrics(res, stages_):
    """Per (organism, vector): projection, cosine, |CAP|, proj_rank, vector_main_phase,
    plus each organism's best-aligned driver. Returns (long_df, best_driver_series)."""
    S = res['species'][['CAP1', 'CAP2']]
    mag = np.hypot(S['CAP1'], S['CAP2'])
    sites = res['sites'][['CAP1', 'CAP2']].copy(); sites['ph'] = stages_
    cent = sites.groupby('ph')[['CAP1', 'CAP2']].mean()
    cos_by_vec, rows = {}, []
    for v, row in res['biplot'][['CAP1', 'CAP2']].iterrows():
        vhat = row.values / (np.linalg.norm(row.values) + 1e-12)
        proj = pd.Series(S.values @ vhat, index=S.index)
        cos = proj / (mag + 1e-12)
        cos_by_vec[v] = cos
        vmain = max({ph: float(cent.loc[ph].values @ vhat) for ph in cent.index}.items(),
                    key=lambda x: x[1])[0]
        rank = proj.rank(ascending=False).astype(int)
        for g in S.index:
            rows.append(dict(organism=g, vector=v, vector_main_phase=vmain,
                             projection=round(float(proj[g]), 4), cosine=round(float(cos[g]), 4),
                             abs_CAP=round(float(mag[g]), 4), proj_rank=int(rank[g])))
    best = pd.DataFrame(cos_by_vec).idxmax(axis=1)
    return pd.DataFrame(rows), best


def run_threshold(suffix, Y, sample_ids, stages, peak, top29):
    panels = build_panels(Y, sample_ids, stages)
    # ---- evidence CSV (both panels) ----
    ev_frames = []
    metrics = {}
    for panel, res, st in panels:
        mod = M.compute_modules(res['species'][['CAP1', 'CAP2']])['louvain']
        df, best = panel_metrics(res, st)
        metrics[panel] = (res, st, df, best)
        df = df.assign(panel=panel,
                       module=df['organism'].map(lambda g: f"M{int(mod.get(g, 0))}"),
                       peak_phase=df['organism'].map(lambda g: peak.get(g, '')),
                       in_prev_top10=df['organism'].map(lambda g: g in top29),
                       is_best_driver=df.apply(lambda r: best.get(r['organism']) == r['vector'], axis=1))
        ev_frames.append(df)
    ev = pd.concat(ev_frames, ignore_index=True)
    ev = ev[['organism', 'module', 'peak_phase', 'in_prev_top10', 'panel', 'vector',
             'vector_main_phase', 'projection', 'cosine', 'abs_CAP', 'proj_rank', 'is_best_driver']]
    ev = ev.sort_values(['panel', 'vector', 'projection'], ascending=[True, True, False])
    ev_path = OUT / f'dbRDA_evidence_genus_vector{suffix}.csv'
    ev.to_csv(ev_path, index=False)
    print(f'wrote {ev_path}  ({ev.organism.nunique()} organisms x {ev.vector.nunique()} vectors, {len(ev)} rows)')

    # ---- association tables ----
    if suffix == '':
        _write_base_assoc()
    else:
        _write_5pct_assoc(suffix, metrics, top29)


def _write_base_assoc():
    cols = ["vector", "strongest_associated_genera", "main_phase", "module", "interpretation"]
    env = [(v, CURATED_ENV_GENERA[v], ph, mod, interp) for (v, ph, mod, interp) in CURATED_ENV]
    perf = [(v, CURATED_PERF_GENERA[v], ph, mod, interp) for (v, ph, mod, interp) in CURATED_PERF]
    perf += [
        ("Opposite high-performance vectors",
         "Thauera; Zoogloea; Sediminibacterium; Kapabacteriales; Pseudoxanthomonas", "I", "M1",
         "The early aerobic community is the geometric inverse of all later denitrification/removal/N2O states (strongly negative projections)."),
        ("Weak / near-origin",
         "Aeromonas; Hydrogenophaga; Ignavibacteriales; Dokdonella; Holophagae", "Transitional", "mixed",
         "Weak/non-specific performance association (smallest |CAP|); Aeromonas's faint signal points to N2O."),
    ]
    pd.DataFrame(env, columns=cols).to_csv(OUT / 'dbRDA_environmental_vector_associations.csv', index=False)
    pd.DataFrame(perf, columns=cols).to_csv(OUT / 'dbRDA_performance_vector_associations.csv', index=False)
    print(f"wrote base environmental ({len(env)}) + performance ({len(perf)}) association tables")


def _top_for_vector(df, vraw, top29, n=TOPN_5PCT, cos_min=COS_MIN, sign=+1):
    """Top genera for one vector by signed projection, cos-filtered; '*' marks
    genera beyond the original top-10 union."""
    d = df[df['vector'] == vraw].copy()
    d['s'] = sign * d['projection']
    d = d[(d['s'] > 0) & (sign * d['cosine'] >= cos_min)].sort_values('s', ascending=False).head(n)
    return "; ".join(disp(g) + ("" if g in top29 else "*") for g in d['organism'])


def _write_5pct_assoc(suffix, metrics, top29):
    cols = ["vector", "strongest_associated_genera_5pct", "main_phase", "module", "interpretation"]
    _, _, dfE, _ = metrics['environment']
    _, _, dfP, _ = metrics['performance']

    # Carry the base main phase / interpretation (these tables expand the originals);
    # the 5%-fit per-vector main phase is recorded in the evidence CSV instead.
    env = []
    for (v, ph_base, mod, interp) in CURATED_ENV:
        genera = _top_for_vector(dfE, ENV_VEC_RAW.get(v, v), top29)
        env.append((v, genera, ph_base, mod, interp))

    perf = []
    for (v, ph_base, mod, interp) in CURATED_PERF:
        genera = _top_for_vector(dfP, PERF_VEC_RAW[v], top29)
        perf.append((v, genera, ph_base, mod, interp))
    # special rows over the broader pool
    # opposite = most negative MEAN projection across the 4 performance vectors
    meanproj = dfP.groupby('organism')['projection'].mean().sort_values()
    opp = "; ".join(disp(g) + ("" if g in top29 else "*") for g in meanproj.head(TOPN_5PCT).index)
    perf.append(("Opposite high-performance vectors", opp, "I", "M1",
                 "The early aerobic community is the geometric inverse of all later denitrification/removal/N2O states (strongly negative projections)."))
    # weak = smallest |CAP| (least specific association)
    magP = dfP.groupby('organism')['abs_CAP'].first().sort_values()
    weak = "; ".join(disp(g) + ("" if g in top29 else "*") for g in magP.head(TOPN_5PCT).index)
    perf.append(("Weak / near-origin", weak, "Transitional", "mixed",
                 "Weak/non-specific performance association (smallest |CAP|)."))

    pd.DataFrame(env, columns=cols).to_csv(OUT / f'dbRDA_environmental_vector_associations{suffix}.csv', index=False)
    pd.DataFrame(perf, columns=cols).to_csv(OUT / f'dbRDA_performance_vector_associations{suffix}.csv', index=False)
    print(f"wrote _5% environmental ({len(env)}) + performance ({len(perf)}) association tables  "
          f"('*' = genus beyond the original top-10 union)")


def main():
    peak = peak_phase_map(M.LONG_FILE)
    # base: top-10-per-phase union (29 organisms)
    Yb, sidb, stb, top10b = M.build_abundance_matrix(M.LONG_FILE)
    top29 = set(Yb.columns)
    print(f"=== base threshold: {Yb.shape[1]} organisms ===")
    run_threshold('', Yb, sidb, stb, peak, top29)
    # 5%: max per-sample rel. abundance >= 5%
    Y5, sid5, st5, _ = M.build_abundance_matrix(M.LONG_FILE, max_rel_threshold=0.05)
    print(f"\n=== 5% threshold: {Y5.shape[1]} organisms ({Y5.shape[1]-len(top29)} beyond the top-10 union) ===")
    run_threshold('_5%', Y5, sid5, st5, peak, top29)


if __name__ == '__main__':
    main()
