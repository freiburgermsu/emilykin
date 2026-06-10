#!/usr/bin/env python
"""
validate_edlib_recall.py — does the edlib top-500 prefilter lose any true biopython hits?

For ~20 ASVs stratified across the best-identity range, run an EXHAUSTIVE biopython local
Smith-Waterman search over EVERY unique BV-BRC 16S reference (same scoring scheme and same
reference universe as the production pipeline), keep the top-5, then measure how many of those
exhaustive top-5 are captured by the pipeline's edlib->biopython top-20.

Because both searches use identical biopython scoring on the identical 459K-ref universe, any
gap is attributable solely to the edlib top-500 candidate prefilter.

Outputs (in --outdir):
  asv_top5_alignment_hits_validation.json   exhaustive top-5 per ASV, same hit schema as
                                            asv_top20_alignment_hits.json (+ in_pipeline_top20,
                                            pipeline_rank annotations)
  validation_recall_per_asv.csv             per-ASV capture stats
  validation_recall_summary.{csv,md}        aggregate statistics (headline % captured)
"""
from __future__ import annotations
import argparse, csv, heapq, json, os, time
from multiprocessing import Pool

import numpy as np
import edlib

import edlib_biopython_hits as P   # reuse loaders, scoring scheme, header/lineage helpers

# ----- fork-shared globals (set in main before Pool) -----
REF_BLOB = b""
REF_OFF = None
ASV_SEQS = []      # the validation ASV sequences (<= ~20)
KEEP = 20          # per-chunk / merge depth (>= 5, headroom for ties)
ALIGNER = None


def init_worker():
    global ALIGNER
    ALIGNER = P.make_aligner()


def score_chunk(span):
    """Score refs[lo:hi] against every validation ASV; return per-ASV local top-KEEP (score, ref_idx)."""
    lo, hi = span
    na = len(ASV_SEQS)
    heaps = [[] for _ in range(na)]      # min-heaps of (score, ref_idx)
    for i in range(lo, hi):
        ref = REF_BLOB[REF_OFF[i]:REF_OFF[i + 1]].decode("ascii")
        for a in range(na):
            sc = ALIGNER.score(ASV_SEQS[a], ref)
            h = heaps[a]
            if len(h) < KEEP:
                heapq.heappush(h, (sc, i))
            elif sc > h[0][0]:
                heapq.heapreplace(h, (sc, i))
    return heaps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 4))
    ap.add_argument("--n-asvs", type=int, default=20)
    ap.add_argument("--topk", type=int, default=5, help="exhaustive hits kept per ASV")
    ap.add_argument("--min-ref-len", type=int, default=200)
    ap.add_argument("--chunk", type=int, default=512)
    ap.add_argument("--pipeline-json", default="asv_top20_alignment_hits.json")
    ap.add_argument("--summary-csv", default="asv_alignment_summary.csv")
    args = ap.parse_args()
    global REF_BLOB, REF_OFF, ASV_SEQS
    wall0 = time.perf_counter()

    # references (identical universe to the pipeline)
    print("[load] references ...", flush=True)
    md5s, REF_BLOB, REF_OFF, ref_len = P.load_references(args.min_ref_len)
    n_ref = len(md5s)
    md5_index = {m: i for i, m in enumerate(md5s)}
    print(f"[load] {n_ref:,} unique refs", flush=True)

    # ASV sequences + the pipeline results to validate against
    asv_ids, asv_seqs = P.load_asvs(P.ASVS_FASTA)
    seq_of = dict(zip(asv_ids, asv_seqs))
    midas = P.load_midas(P.TAXONOMY_CSV)
    pipe = json.load(open(os.path.join(args.outdir, args.pipeline_json)))

    # stratified pick: sort by pipeline best_identity, take n evenly spaced by rank
    rows = list(csv.DictReader(open(os.path.join(args.outdir, args.summary_csv))))
    rows = [r for r in rows if r["best_identity"]]
    rows.sort(key=lambda r: float(r["best_identity"]))
    N = len(rows)
    n = min(args.n_asvs, N)
    pick_idx = [round(i * (N - 1) / (n - 1)) for i in range(n)]
    picked = [rows[i]["asv"] for i in pick_idx]
    ASV_SEQS = [seq_of[a] for a in picked]
    print(f"[pick] {n} ASVs stratified over best_identity "
          f"[{rows[0]['best_identity']} .. {rows[-1]['best_identity']}]", flush=True)

    # ---- exhaustive biopython search, parallel over ref chunks ----
    spans = [(lo, min(lo + args.chunk, n_ref)) for lo in range(0, n_ref, args.chunk)]
    merged = [[] for _ in range(n)]      # per-ASV accumulated (score, ref_idx)
    t0 = time.perf_counter()
    with Pool(args.workers, initializer=init_worker) as pool:
        done = 0
        for heaps in pool.imap_unordered(score_chunk, spans, chunksize=1):
            for a in range(n):
                merged[a].extend(heaps[a])
            done += 1
            if done % 100 == 0 or done == len(spans):
                el = time.perf_counter() - t0
                print(f"[scan] chunk {done}/{len(spans)}  elapsed={el:5.0f}s  "
                      f"eta={el/done*(len(spans)-done):5.0f}s", flush=True)
    scan_wall = time.perf_counter() - t0

    # global top-k per ASV (score desc, ref_idx asc as deterministic tiebreak)
    exh = []
    for a in range(n):
        best = sorted(set(merged[a]), key=lambda x: (-x[0], x[1]))[:args.topk]
        exh.append(best)

    # ---- enrichment + capture annotation ----
    print("[enrich] md5->header + taxdump ...", flush=True)
    md5_hdr = json.load(open(P.DB_MD5_ID))
    import taxopy
    taxdb = taxopy.TaxDb(nodes_dmp=P.TAXDUMP_NODES, names_dmp=P.TAXDUMP_NAMES)
    lin_cache = {}

    def lineage_for(tid):
        if tid is None:
            return {r: None for r in P.RANKS}
        if tid not in lin_cache:
            try:
                rd = taxopy.Taxon(tid, taxdb).rank_name_dictionary
                lin_cache[tid] = {"Kingdom": rd.get("superkingdom") or rd.get("kingdom") or rd.get("domain"),
                                  "Phylum": rd.get("phylum"), "Class": rd.get("class"), "Order": rd.get("order"),
                                  "Family": rd.get("family"), "Genus": rd.get("genus"), "Species": rd.get("species")}
            except Exception:
                lin_cache[tid] = {r: None for r in P.RANKS}
        return lin_cache[tid]

    aligner = P.make_aligner()
    validation, per_asv = {}, []
    tot_capt = tot_genuine_miss = full_capture = top1_in20 = top1_in5 = top1_md5_match = top1_score_match = 0

    for a, asv in enumerate(picked):
        aseq = ASV_SEQS[a]
        p = pipe.get(asv, {})
        p_top = p.get("top20", [])
        p_md5_rank = {h["md5"]: h["rank"] for h in p_top}
        p_scores = [h["align_score"] for h in p_top]
        p_min20 = min(p_scores) if p_scores else float("-inf")
        p_top1_md5 = p_top[0]["md5"] if p_top else None
        p_top1_score = p_top[0]["align_score"] if p_top else None

        hits, e5_md5 = [], []
        for rank, (sc, idx) in enumerate(exh[a], 1):
            md5 = md5s[idx]
            e5_md5.append(md5)
            ref = REF_BLOB[REF_OFF[idx]:REF_OFF[idx + 1]].decode("ascii")
            aln = aligner.align(aseq, ref)[0]
            ident, matches, alen, rs, re_ = P._identity_and_coords(aln)
            ed = edlib.align(aseq, ref, mode="HW", task="distance")["editDistance"]
            org, gid, tid, feat = P.parse_header(md5_hdr.get(md5, ""))
            hits.append({
                "rank": rank, "align_score": float(sc), "identity": round(ident, 4),
                "n_matches": matches, "aligned_len": alen,
                "edlib_distance": int(ed), "edlib_identity": round(1 - ed / max(1, len(aseq)), 4),
                "organism": org, "genome_id": gid, "taxon_id": tid, "feature_id": feat, "md5": md5,
                "ref_seq_len": int(ref_len[idx]), "ref_aln_start": rs, "ref_aln_end": re_,
                "lineage": lineage_for(tid),
                "in_pipeline_top20": md5 in p_md5_rank,
                "pipeline_rank": p_md5_rank.get(md5),
            })
        midas_tax, rel_ab = midas.get(asv, ("", 0.0))
        validation[asv] = {
            "asv_len": len(aseq), "midas_taxonomy": midas_tax, "rel_ab": rel_ab,
            "search": "exhaustive_biopython_local_all_BVBRC", "n_refs_searched": n_ref,
            "best_align_score": hits[0]["align_score"] if hits else None,
            "best_identity": hits[0]["identity"] if hits else None,
            "top5": hits,
        }

        # capture stats
        captured = sum(1 for m in e5_md5 if m in p_md5_rank)
        genuine_miss = sum(1 for (sc, idx) in exh[a]
                           if md5s[idx] not in p_md5_rank and sc > p_min20)
        e_top1_md5, e_top1_score = e5_md5[0], exh[a][0][0]
        in20 = e_top1_md5 in p_md5_rank
        in5 = p_md5_rank.get(e_top1_md5, 99) <= 5
        md5_match = (e_top1_md5 == p_top1_md5)
        score_match = (p_top1_score is not None and abs(e_top1_score - p_top1_score) < 1e-6)
        tot_capt += captured; tot_genuine_miss += genuine_miss
        full_capture += (captured == len(e5_md5))
        top1_in20 += in20; top1_in5 += in5; top1_md5_match += md5_match; top1_score_match += score_match
        per_asv.append({
            "asv": asv, "rel_ab": rel_ab, "pipe_best_identity": p.get("best_identity"),
            "midas_genus": midas_tax.split()[-2] if len(midas_tax.split()) >= 2 else "",
            "captured_of_5": captured, "recall_5in20": round(captured / len(e5_md5), 3),
            "genuine_misses": genuine_miss,
            "exh_top1_score": round(e_top1_score, 1), "pipe_top1_score": p_top1_score,
            "top1_score_match": score_match, "top1_md5_match": md5_match,
            "exh_top1_in_pipe_top20": in20, "exh_top1_in_pipe_top5": in5,
        })

    # ---- write outputs ----
    json.dump(validation, open(os.path.join(args.outdir, "asv_top5_alignment_hits_validation.json"), "w"), indent=1)
    with open(os.path.join(args.outdir, "validation_recall_per_asv.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_asv[0].keys())); w.writeheader(); w.writerows(per_asv)

    total_e5 = sum(len(e) for e in exh)
    summary = {
        "n_asvs_validated": n,
        "exhaustive_top5_hits_total": total_e5,
        "pct_top5_captured_in_pipeline_top20": round(100 * tot_capt / total_e5, 2),
        "n_genuine_misses_score_better_than_pipeline_20th": tot_genuine_miss,
        "pct_asvs_full_5of5_capture": round(100 * full_capture / n, 2),
        "pct_asvs_exhaustive_top1_in_pipeline_top20": round(100 * top1_in20 / n, 2),
        "pct_asvs_exhaustive_top1_in_pipeline_top5": round(100 * top1_in5 / n, 2),
        "pct_asvs_exhaustive_top1_is_pipeline_top1_md5": round(100 * top1_md5_match / n, 2),
        "pct_asvs_top1_score_matches_pipeline": round(100 * top1_score_match / n, 2),
        "exhaustive_search_wall_seconds": round(scan_wall, 1),
        "total_wall_seconds": round(time.perf_counter() - wall0, 1),
    }
    json.dump(summary, open(os.path.join(args.outdir, "validation_recall_summary.csv").replace(".csv", ".json"), "w"), indent=2)
    # markdown + csv table of the summary
    with open(os.path.join(args.outdir, "validation_recall_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["metric", "value"])
        for k, v in summary.items():
            w.writerow([k, v])
    with open(os.path.join(args.outdir, "validation_recall_summary.md"), "w") as fh:
        fh.write("# edlib-prefilter recall validation\n\n")
        fh.write("Exhaustive biopython local search over all unique BV-BRC 16S refs vs. the "
                 "edlib(top-500)->biopython(top-20) pipeline, identical scoring and reference universe.\n\n")
        fh.write("| metric | value |\n|---|---|\n")
        for k, v in summary.items():
            fh.write(f"| {k} | {v} |\n")
    print("\n=== VALIDATION SUMMARY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
