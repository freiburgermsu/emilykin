# BV-BRC embedding hits for the EmilyKin EBPR ASVs

This folder maps the study's 16S **V4–V5 ASVs** (`../qiime/`) against a BV-BRC 16S **embedding space**
built with Nucleotide-Transformer-v2-500m — the *same* embedding mapping and report/JSON generation
used for the anaerobic-digester study (`~/Documents/prFBA`). The reference space is **region-matched**:
BV-BRC 16S genes had their V4–V5 insert excised by in-silico PCR (515F/926R) and embedded identically,
so amplicons compare like-with-like. **Read `REPORT_emilykin_hits.md` for the findings.**

## Headline

The 3,950 ASVs (77 samples) are placed confidently by sequence — **median best cosine 0.9985**,
98.7% ≥ 0.97, 791 at exact 1.0 — but **taxonomic-name transfer is weak** (abundance-weighted genus
concordance 0.40, family **0.16**), the *inverse* of the digester study (~0.89). The reason is
ecological: EmilyKin's dominant members are **uncultured candidate PAOs/GAOs** (Ca_Competibacter,
Ca_Accumulibacter, …) whose near-exact BV-BRC matches are **unclassified MAGs with no family name**.
Sequence placement is reliable for both communities; *name* transfer depends on whether the abundant
taxa are cultured. Any-top-20 concordance (genus 0.66) and the cosine itself are the fair signals.

## Contents

| | |
|---|---|
| `REPORT_emilykin_hits.md` | findings + adversarial-verification appendix |
| **`asv_top20_hits.json`** | **per-ASV top-20 BV-BRC hits**: cosine + organism/genome_name/taxon_id/genome_id/feature_id/n_genomes |
| `asv_summary.csv` | per-ASV best hit, cosine, naive cosine, concordance, novelty |
| `asv_concordance.csv` / `concordance_by_rank.json` | multi-rank concordance vs MiDAS (Ca_/synonym-normalized) |
| `findings_stats.json` | aggregate match statistics |
| `prep_inputs.py` | extract qiime artifacts → `inputs/` (asvs.fasta, taxonomy.csv, asv_IDs.csv, per-ASV rel_ab) |
| `hit_amplicons.py` | embed ASVs → top-20 cosine hits + BV-BRC metadata (+ naive comparison) |
| `concordance.py` | taxonomic concordance vs MiDAS (NCBI taxdump via taxopy) |
| `nt_embed.py` | shared NT-v2 encoder |
| `inputs/` | derived inputs (small text files; `feature-table.biom` git-ignored, regenerable) |

## Reproduce

```bash
# from this directory; venv = ~/Documents/py_venv
python prep_inputs.py                                   # ../qiime/*.qza -> inputs/
PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True python hit_amplicons.py \
    --amplicons inputs/asvs.fasta \
    --store      ~/Documents/prFBA/v4v5_store \          # region-matched BV-BRC V4-V5 refs (reused)
    --naive-store ~/Documents/prFBA \                    # full-length refs (naive comparison)
    --asv-ids inputs/asv_IDs.csv --taxonomy inputs/taxonomy.csv --topk 20 --outdir .
python concordance.py --taxonomy inputs/taxonomy.csv --hits asv_top20_hits.json --outdir .
```

The region-matched reference store (`prFBA/v4v5_store`) and the NCBI taxdump are reused from the
prFBA project (large, git-ignored there); rebuild them with prFBA's `insilico_pcr.py` + `embed_16s.py`
if absent.
