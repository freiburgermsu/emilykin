"""Top-20 overlap: edlib-prefilter CPU pipeline vs prefilter-free GPU exhaustive.

Overlap is the per-ASV Jaccard-style set agreement of the two top-20 hit sets (by md5):
    overlap% = |CPU_top20 ∩ GPU_top20| / 20.
Also reports a *score-aware* overlap that forgives pure tie-swaps (a non-shared hit whose
SW score equals one already in the CPU top-20), to separate genuine recall loss from
boundary tie-breaking noise. Both top-20s are ranked by the identical local-SW score.
"""
import json, statistics as st
from collections import Counter

D = "/home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/"
cpu = json.load(open(D + "asv_top20_alignment_hits.json"))
gpu = json.load(open(D + "asv_top500_alignment_hits_gpu.json"))

rows = []
for asv, ce in cpu.items():
    ge = gpu.get(asv)
    if not ge:
        continue
    c = ce["top20"]
    g = ge["top500"][:20]                                  # GPU top-20 (rank-ordered)
    cset = {h["md5"] for h in c}
    gset = {h["md5"] for h in g}
    inter = cset & gset
    raw = len(inter) / 20.0
    # score-aware: count a GPU hit as covered if CPU top-20 holds the same score (tie-swap)
    cscores = Counter(round(h["align_score"], 3) for h in c)
    covered = 0
    for h in g:
        if h["md5"] in cset:
            covered += 1
        elif cscores.get(round(h["align_score"], 3), 0) > 0:
            covered += 1
            cscores[round(h["align_score"], 3)] -= 1          # consume the tie slot
    saware = covered / 20.0
    rows.append({
        "asv": asv, "best_id": ge.get("best_identity") or 0.0, "rel_ab": ce.get("rel_ab", 0.0),
        "raw": raw, "saware": saware, "inter": len(inter),
        "midas": ce.get("midas_taxonomy", "").split(" ")[-1] if ce.get("midas_taxonomy") else "",
    })

n = len(rows)
raws = [r["raw"] for r in rows]
saws = [r["saware"] for r in rows]
print(f"ASVs compared: {n}\n")
print("================  TOP-20 OVERLAP: CPU prefilter pipeline vs GPU exhaustive  ================")
print(f"  mean raw set overlap          : {100*st.mean(raws):.3f}%   (= {st.mean(raws)*20:.2f} of 20 hits shared on average)")
print(f"  median raw set overlap        : {100*st.median(raws):.3f}%")
print(f"  mean score-aware overlap      : {100*st.mean(saws):.3f}%   (tie-swaps forgiven)")
print(f"  ASVs with 20/20 identical     : {sum(1 for r in rows if r['inter']==20)}  ({100*sum(1 for r in rows if r['inter']==20)/n:.2f}%)")
print(f"  ASVs with >=19/20             : {sum(1 for r in rows if r['inter']>=19)}  ({100*sum(1 for r in rows if r['inter']>=19)/n:.2f}%)")
print(f"  ASVs with >=15/20             : {sum(1 for r in rows if r['inter']>=15)}  ({100*sum(1 for r in rows if r['inter']>=15)/n:.2f}%)")
print(f"  ASVs with < 10/20             : {sum(1 for r in rows if r['inter']<10)}  ({100*sum(1 for r in rows if r['inter']<10)/n:.2f}%)")

print("\n---- distribution of shared hits (out of 20) ----")
hist = Counter(r["inter"] for r in rows)
for k in range(20, -1, -1):
    if hist.get(k):
        bar = "#" * max(1, round(60 * hist[k] / n))
        print(f"  {k:2d}/20 : {hist[k]:4d} ASVs ({100*hist[k]/n:5.2f}%) {bar}")

print("\n---- mean overlap stratified by ASV best identity (difficulty) ----")
bins = [(0,0.80,"<0.80"),(0.80,0.90,"0.80-0.90"),(0.90,0.95,"0.90-0.95"),(0.95,0.97,"0.95-0.97"),
        (0.97,0.99,"0.97-0.99"),(0.99,0.999,"0.99-0.999"),(0.999,1.01,">=0.999")]
print(f"  {'identity':>11} | {'n':>5} | {'mean raw%':>9} | {'mean score-aware%':>17} | {'%ASVs 20/20':>11}")
for lo, hi, lab in bins:
    sub = [r for r in rows if lo <= r["best_id"] < hi]
    if not sub:
        continue
    full = sum(1 for r in sub if r["inter"] == 20)
    print(f"  {lab:>11} | {len(sub):>5} | {100*st.mean(r['raw'] for r in sub):>8.2f}% | "
          f"{100*st.mean(r['saware'] for r in sub):>16.2f}% | {100*full/len(sub):>10.1f}%")

print("\n---- 12 worst-overlap ASVs ----")
for r in sorted(rows, key=lambda x: (x["raw"], -x["best_id"]))[:12]:
    print(f"  {r['asv'][:12]} shared={r['inter']:2d}/20 raw={100*r['raw']:5.1f}% score-aware={100*r['saware']:5.1f}% "
          f"best_id={r['best_id']:.3f} relab={r['rel_ab']:.4f} midas={r['midas']}")

# save concise artifacts
with open(D + "top20_overlap_per_asv.csv", "w") as fh:
    fh.write("asv,best_identity,rel_ab,shared_of_20,raw_overlap_pct,score_aware_overlap_pct,midas_tail\n")
    for r in sorted(rows, key=lambda x: x["raw"]):
        fh.write(f"{r['asv']},{r['best_id']:.4f},{r['rel_ab']:.6f},{r['inter']},"
                 f"{100*r['raw']:.2f},{100*r['saware']:.2f},{r['midas']}\n")
print(f"\n[wrote] top20_overlap_per_asv.csv ({n} rows)")
