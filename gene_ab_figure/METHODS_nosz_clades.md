# nosZ (K00376) Clade I / Clade II classification — methods

Implemented in `05_classify_nosz.py`; outputs `data/nosz_clades.tsv` (per-gene
scores + confidence) and feeds the split nosZ rows of the figure.

## 1. Reference database (the source materials in `literature_claudeI_II/`)

Nucleotide nosZ reference **alignments** (aligned FASTA), used as the two-clade
reference set:

| File | Clade | n seqs | alignment columns | representative taxa |
|------|-------|-------:|------------------:|---------------------|
| `NosZREF-CladeI`  | **I** (typical / denitrifier) | 157 | 254 | *Pseudomonas, Bradyrhizobium, Rhodanobacter, Neisseria, Azoarcus, Thauera, Achromobacter…* (broad Alpha/Beta/Gammaproteobacteria) |
| `NosZREF-1577A` | II.A | 80 | 373 | Gemmatimonadetes + environmental clones |
| `NosZREF-1577B` | II.B | 72 | 378 | *Prevotella, Gemmatimonas* |
| `NosZREF-1577C` | II.C | 71 | 374 | Bacteroidetes (*Prevotella, Dyadobacter, Runella, Haliscomenobacter*) |
| `NosZREF-1577D` | II.D | 34 | 546 | *Dechloromonas, Dechlorosoma* |
| `NosZREF_1577E` | II.E | 60 | 375 | *Desulfomonile, Desulfitobacterium, Salinibacter* |
| `NosZREF_1577G` | II.G | 34 | 526 | *Desulfitobacterium* + environmental clones |
| `NosZREF-1577H` | II.H | 23 | 410 | *Desulfosporosinus, Desulfitobacterium* (Firmicutes) |

Clade II total = **374** sequences across 7 sub-clades; reference total = **531**.

### Modifications made to the source materials

- **Reference sequence content — NONE.** The aligned FASTAs were consumed
  verbatim. No realignment, degapping, trimming, translation, dereplication, or
  subsetting; every sequence in every file was used.
- **Clade grouping — an interpretation, not an edit.** The single `CladeI` file
  was taken as Clade I; the seven `1577A–H` files were treated as Clade II
  sub-clades. The folder carried no explicit I/II key — this grouping was
  inferred from (a) the file naming (`CladeI` vs `1577*`) and (b) the taxon
  composition of the `1577*` sets, which are the canonical atypical /
  non-denitrifier Clade II lineages (Bacteroidetes, Gemmatimonadetes,
  Desulfo* Firmicutes, *Dechloromonas*, *Salinibacter*).
- **Folder housekeeping (non-content).** The folder was renamed from
  `literature_claudeI␣II` (the separator was a stray private-use Unicode
  character, U+F027) to `literature_claudeI_II`, and macOS AppleDouble sidecar
  files (`._*`) were deleted. No `.fasta.text` payload was touched.

## 2. Query sequences (our MAGs)

- The 13 nosZ genes are every K00376 locus in `data/target_genes.bed` (one per
  MAG).
- Each CDS was extracted from `data/combined_ref.fa` by slicing the contig at the
  BED `[start, end)` interval (0-based, half-open), upper-cased (U→T). Strand was
  **not** resolved — `nhmmer` searches both strands, so it does not affect the
  call.

## 3. Assignment software and procedure

- **Aligner / search tool: HMMER3 `nhmmer`** (nucleotide profile homology
  search), run through **pyhmmer v0.12.1** (Python bindings — no standalone
  HMMER binary was installed; no BLAST/DIAMOND/MMseqs were used). `pyrodigal`
  was available but not used.
- One profile **HMM was built per reference alignment** with
  `pyhmmer.plan7.Builder` under default HMMER3 model-construction settings
  (DNA alphabet, default background; symfrac-based match-state assignment, so a
  few highly-gapped columns become inserts → HMM lengths 248/364/373/371/486/
  367/480/378; Henikoff position-based sequence weighting).
- `nhmmer` scored each of the 13 genes against all 8 clade HMMs (both strands);
  the best bit-score per clade was kept.
- **Decision rule:** Clade I if (Clade-I bit-score) ≥ max(seven Clade-II
  sub-clade bit-scores), else Clade II.
- **Confidence:** `high` if |score margin| ≥ 40 bits, else `low` (flagged with
  ○ on the figure and `confidence=low` in `data/nosz_clades.tsv`).

## 4. Result

12 / 13 = **Clade II** (10 unambiguous, margins −195 to −366 bits;
`CAN_5_bin.40` *Giesbergeria* low-confidence at −24). 1 = **Clade I**,
`coasm_bin.55` *Competibacter*, low-confidence at **+7.9 bits** (a near-tie:
195.4 vs 187.5). Full per-gene scores in `data/nosz_clades.tsv`.

## 5. Distinction from the abundance (RPKM) step

The clade assignment above operates on **assembled gene sequences** with
`nhmmer`. The **RPKM abundances** shown in the figure come from a separate,
earlier step: long Nanopore reads were mapped to the dereplicated-MAG reference
with **minimap2 (`-ax map-ont`, `--secondary=no`)** on the cluster
(`02_build_and_align.sh`) and counted per gene. minimap2 was the read→MAG
aligner; `nhmmer` was the gene→clade classifier. They are independent.
