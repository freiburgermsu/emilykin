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
