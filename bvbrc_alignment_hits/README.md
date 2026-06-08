# BV-BRC alignment hits for the EmilyKin EBPR ASVs (edlib → biopython)

This folder maps the study's 16S **V4–V5 ASVs** (`../asvs.fasta`, 3,950 sequences) against the
BV-BRC 16S reference database by **sequence alignment**, as a complement to the embedding-based
placement in `../bvbrc_embedding_hits/`. Two stages, fully parallelized across all cores:

1. **edlib (infix / `HW` mode), all references.** For each ASV, the edit distance of its best
   placement *as an infix* within every unique BV-BRC 16S gene is computed, and the **top 500**
   lowest-distance references are kept. An adaptive `k`-band (edlib `k` = the running 500-th best
   distance) lets the >99% of clearly-dissimilar references reject in microseconds, so the full
   3,950 × ~459K scan is tractable.
2. **biopython `PairwiseAligner` local Smith-Waterman, 500 candidates.** Those 500 candidates are
   rescored with a BLASTN-like nucleotide scheme; the **top 20** by alignment score are kept and
   re-aligned to record % identity, aligned length, and reference coordinates.

The reference set is the **deduplicated** unique-sequence set (`16S_md5_seq.json`, md5 → sequence;
~459K unique 16S genes after dropping fragments < 200 bp). Identical reference sequences collapse to
a single hit; each is mapped back to a representative BV-BRC genome via `16S_md5_ID.json`, with
`taxon_id` parsed from the genome_id and a full NCBI lineage attached via `taxopy` (prFBA taxdump).

## Parameters

| stage | tool | setting |
|---|---|---|
| 1 | edlib | `mode="HW"` (query=ASV as infix of target=gene), `task="distance"`, adaptive `k`, top **500** |
| 2 | biopython | `PairwiseAligner` `mode="local"`; match **+2**, mismatch **−3**, gap_open **−5**, gap_extend **−2**; top **20** |

## Outputs

| file | contents |
|---|---|
| **`asv_top20_alignment_hits.json`** | **the mapping JSON** — `{asv_hash: {asv_len, midas_taxonomy, rel_ab, best_align_score, best_identity, n_edlib_candidates, top20:[…]}}`. Each `top20` hit: `rank, align_score, identity, n_matches, aligned_len, edlib_distance, edlib_identity, organism, genome_id, taxon_id, feature_id, md5, ref_seq_len, ref_aln_start, ref_aln_end, lineage{Kingdom…Species}` |
| `asv_alignment_summary.csv` | one row per ASV: best-hit organism/genome/score/identity + best genus/family |
| `run_stats.json` | parameters, counts, **per-stage CPU seconds** (edlib vs biopython) and **wall time** |
| `edlib_biopython_hits.py` | the pipeline (this is what produced the above) |
| `run.log` | progress log of the production run |

`identity` = matching columns / aligned columns from the biopython local alignment;
`edlib_identity` = 1 − edit_distance / asv_len (the stage-1 score the candidate was ranked by).

## Reproduce

```bash
# venv = ~/Documents/py_venv ; run from this directory
python edlib_biopython_hits.py --workers 60 --chunksize 1 --k1 500 --k2 20 --outdir .
# --limit N   benchmarks on the first N ASVs (the full reference load still happens)
```

Inputs are read directly from `~/Documents/codiffusion_bioreactor/model_inputs/`
(`16S_md5_seq.json`, `16S_md5_ID.json`) and the NCBI taxdump at `~/Documents/prFBA/{nodes,names}.dmp`;
MiDAS lineage + per-ASV max relative abundance are reused from `../bvbrc_embedding_hits/inputs/taxonomy.csv`.

## Prefilter validation

`validate_edlib_recall.py` tests whether the edlib top-500 prefilter loses any true Biopython hit:
for 20 ASVs stratified across the best-identity range it runs an **exhaustive** Biopython local search
over all 459 K refs (top-5) and checks capture in the pipeline top-20.

| file | contents |
|---|---|
| `asv_top5_alignment_hits_validation.json` | exhaustive top-5 per ASV (same hit schema; each hit annotated `in_pipeline_top20`/`pipeline_rank`) |
| `validation_recall_{per_asv.csv,summary.md/csv/json}` | capture statistics |

**Result: 100% recall** — all 100 exhaustive top-5 hits captured in the pipeline top-20, 0 genuine
misses, 20/20 top-1 score matches, even for the hardest ASV at 67.5% identity.

## GPU tier (optional, fastest)

> **The GPU acceleration code now lives in the [`prFBA`](https://github.com/freiburgermsu/prFBA)
> repo** — `gpu_align.py` (kernel + driver), `edlib_biopython_hits.py` (its CPU-pipeline dependency,
> copied so it imports standalone there), and `GPU_ALIGNMENT_FINDINGS.md`. The run outputs it generated
> on **this** data stay in this folder.

`gpu_align.py` — a custom local Smith-Waterman CUDA kernel (CuPy/NVRTC, **no nvcc**) with the **same
scoring scheme** as the Biopython stage, **verified bit-exact** (400/400 vs Biopython) and able to run
**prefilter-free exhaustive** scoring over all 459 K refs. On the RTX 5070 Ti it does the 20-ASV
exhaustive in 66 s @ 57 GCUPS (**7.8× faster than CPU**, 100% reproducing the CPU result); full
3,950-ASV exhaustive extrapolates to ~3.6 h GPU vs ~28 h CPU. See **`GPU_ALIGNMENT_FINDINGS.md`** (in prFBA) for
the full evaluation (incl. why MMseqs2-GPU and CUDASW++4.0 — both **protein-only on GPU** — do not fit
this nucleotide task). Outputs: `gpu_align_stats.json`.

## Full-scale GPU run, comparison, and files not committed to git

The GPU kernel was run **prefilter-free over all 3,950 ASVs × 459,301 refs** (top-500 each,
3.6 h @ 57 GCUPS) and compared against the CPU edlib-prefilter pipeline:

| file | contents |
|---|---|
| `gpu_exhaustive_full_stats.json` | full-run parameters + CPU-vs-pipeline agreement summary |
| `gpu_vs_pipeline_full_comparison_per_asv.csv` | per-ASV GPU-vs-pipeline agreement |
| `compare_top20_overlap.py` → `top20_overlap_per_asv.csv` | per-ASV top-20 set overlap (raw + score-aware) |

**Top-20 overlap (CPU prefilter vs GPU exhaustive): 97.6% mean, 90.9% of ASVs identical, top-1 ~99%.**
The edlib top-500 prefilter is lossy *at depth* for **~5.5% of ASVs** (mostly low-identity / eukaryotic /
off-target), where edit-distance ranking diverges from local-SW score — so the **"100% recall" above is
the 20-ASV validation sample, not the full set**. Use the GPU exhaustive output when deep ranks matter.

> **Not committed to git** (exceeds GitHub's 100 MB file limit; regenerable):
> `asv_top500_alignment_hits_gpu.json` — **1.16 GB**, the full prefilter-free GPU exhaustive top-500
> per ASV (same per-hit schema as `asv_top20_alignment_hits.json`). Regenerate with `python gpu_align.py`
> (now in the prFBA repo; ~3.6 h on the RTX 5070 Ti). It is `.gitignore`d; all other outputs in this folder are tracked.

## Large files (> 100 MB, git-ignored)

| file | size | what it is | regenerate |
|---|---|---|---|
| `asv_top500_alignment_hits_gpu.json` | 1.1 GB | GPU exhaustive Smith-Waterman top-500 hits per ASV | `python gpu_align.py` (in prFBA, exhaustive mode) |

Exceeds GitHub's 100 MB limit → git-ignored (regenerable). Threshold = 100 MB (GitHub hard limit); no file is near 100 GB.
