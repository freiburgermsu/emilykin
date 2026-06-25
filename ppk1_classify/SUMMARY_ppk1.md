# ppk1-based classification of CAN 'Ca. Accumulibacter' MAGs

Mirrors the Petriglieri et al. 2022 (mBio) reevaluation / McDaniel et al. ppk1_Database
workflow: identify ppk1 by HMM, build a nucleotide ppk1 tree with the curated reference
set + outgroup, and read each MAG's ppk1 clade off the tree; corroborate by % identity.

## Inputs & method
- **Reference DB:** elizabethmcd/ppk1_Database (cloned) — `ppk1.hmm`, 40 genome + 41 clone
  reference ppk1 sequences with clade designations (`ppk1-database-info.csv`), and the
  bundled outgroup.
- **Outgroup:** `outgroup-ppk1-coding-regions.fasta` = *Dechloromonas aromatica* +
  *Rhodocyclus tenuis* (the McMahon-lab outgroup shipped with the DB). NB: Petriglieri 2022
  itself rooted on three *Azonexus* (formerly *Dechloromonas*) isolates — the same lineage,
  so the root is equivalent.
- **Our ppk1:** 7 GTDB-*Accumulibacter* MAGs screened. Genes called with **pyrodigal**
  (meta), ppk1 found by **pyhmmer hmmsearch** vs `ppk1.hmm` (top hit, E≤1e-50). 5 MAGs carry
  a detectable ppk1; the matching nucleotide ORF was extracted.
    - `CAN_1_bin.98` had ppk1 (E=3e-153) **despite KO matrix K00937=0** (KO false negative).
    - `CAN_4_bin.225` (KO K00937=1) and `coasm_bin.476` had **no** ppk1 hit even at E≤10
      (fragmented/incomplete MAGs) and were excluded.
    - The "2-copy" MAGs (185/250/347) resolve to exactly one true ppk1 (HMM rejects the paralog).
- **Alignment:** MAFFT v7.526 `--auto` (88 nucleotide seqs × 2315 cols).
- **Tree:** IQ-TREE 2.4.0, ModelFinder `-m MFP` → **TVM+F+I+G4**, 1000 UFBoot + 1000 SH-aLRT,
  rooted on the outgroup. (Paper used MAFFT --auto + IQ-TREE MFP → GTR+F+I+G4, 100 bootstraps —
  equivalent.)
- **Corroboration:** pairwise nucleotide % identity (the DB/ANI criterion: ~90–100% within
  clade, ~80% between clades).

## Results — ppk1 clade and proposed species (Petriglieri 2022 Table 1)

| MAG | ppk1 type | proposed species | SH-aLRT/UFBoot | top %id (clade) |
|-----|-----------|------------------|:--------------:|:---------------:|
| coasm_bin.185 | **IA**  | *Ca.* Accumulibacter **regalis**     | 99.5/100 | 100.0% (IA) |
| coasm_bin.250 | **IC**  | *Ca.* Accumulibacter **delftensis**  | 100/100  | 99.8% (IC) |
| coasm_bin.347 | **IIA** | *Ca.* Accumulibacter **phosphatis**  | 81.9/99  | 99.4% (IIA) |
| CAN_4_bin.64  | Type I (unclassified) | — (skipped, unclear) | 100/100 | 99.1% (≈UBA11064, un) |
| CAN_1_bin.98  | Type I (novel/divergent) | — (skipped, unclear) | — | ~88% (no clade) |
| CAN_4_bin.225 | no ppk1 detected | — | — | — |
| coasm_bin.476 | no ppk1 detected | — | — | — |

**Notes**
- IA→regalis, IIB→propinquus, IIG→affinis are 1:1; IC, IIA (and IIC/IID/IIF) hold multiple
  species, so the species was set from the nearest reference strain (delftensis vs meliphilus
  for IC; phosphatis/UW1 vs aalborgensis for IIA). coasm_bin.250 is 99.8% to *delftensis* and
  coasm_bin.347 is 99.4% to the UW1 *phosphatis* type strain.
- `CAN_1_bin.98` and `CAN_4_bin.64` are Type I but fall outside every named sub-clade (they
  group with the DB's own unclassified *Accumulibacter* sp. UBA11064); species skipped per
  the rule "skip if unclear." `CAN_1_bin.98` is the most divergent (≤88% to any reference) —
  a candidate novel Type I lineage.

## Files
- `ppk1_Database/` — cloned reference database
- `our_ppk1.ffn`, `combined_ppk1.ffn`, `ppk1_aligned.afa` — query/combined seqs + MAFFT alignment
- `ppk1_tree.treefile` / `.iqtree` / `.contree` — IQ-TREE ML tree + report
- `ppk1_tree_annotated.nwk` — tree with clade-labelled tips (open in FigTree/iTOL)
- `ppk1_reference_tree.png/.pdf` — rendered tree, queries highlighted
- `ppk1_classification.tsv` — this table; `partG_placements.json` — placement detail
- `extract_ppk1.py`, `assemble_ppk1.py`, `partG..partI*.py` — pipeline scripts
