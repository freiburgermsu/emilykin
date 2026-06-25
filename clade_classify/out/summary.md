# nosZ Clade Classification — Final Summary

## qNOR / norZ (K11188)
**Completely absent** from all 21 selected MAGs and the entire 337-MAG dataset (zero
hits in KOfamScan or EggNOG). This system exclusively uses the cNOR (norBC) pathway
for NO reduction; no quinol-dependent NOR is present.

## norB (K04561) — correction to figure
norB (the large catalytic subunit of cNOR, paired with norC/K02305) was present in
10/21 selected MAGs above KOfamScan threshold (scores 477–1043) but had been omitted
from the original figure. All 10 norB contigs were already in the BAM reference;
reads were counted from existing alignments. The norC figure row now represents
max(K02305, K04561) RPKM — cNOR coverage increased from 4 to 10 MAGs.

## Part A — K00376 nosZ clade I / II (RESOLVED by tree placement, Part D)

The earlier HMM-vs-diamond disagreement (see `../CLADE_I_II_DISCREPANCY.md`) was
settled phylogenetically: the 13 query proteins were placed in the reference ML
tree (`nosz_tree.nwk`; 330 Chee+Orellana C-NosZ + 259 He refs) and assigned to the
Clade I / Clade II **monophyletic group** they nest in. Clade membership of the
references was set from textbook-unambiguous marker genera only — the disputed
Rhodocyclaceae were deliberately *not* used as anchors, so the tree placed them
independently. Script: `../partD_tree_placement.py`; table: `partD_tree_clades.tsv`.

**CuZ motif:** all 13 = CuA + DXHH → C-NosZ (not L-NosZ / clade III). Confirmed.

**Reconciled clade calls (tree = authority): Clade I = 2, Clade II = 11.**

| MAG | Genus | Clade | nearest reference (sister) | SH support |
|-----|-------|:-----:|----------------------------|:----------:|
| coasm_bin.55 | Competibacter | **I**  | Burkholderia / Rhodoferax / Aromatoleum | 0.97 |
| CAN_5_bin.40 | Giesbergeria  | **I**  | Ralstonia / Cupriavidus / Rubrivivax | 0.67 (borderline) |
| CAN_1_bin.98, coasm_bin.185 | Accumulibacter | II | Dechlorosoma / Dechloromonas | 0.52 |
| CAN_2_bin.6 | Desulfobacillus | II | Dechlorosoma / Dechloromonas | 1.00 |
| CAN_5_bin.70 | Azonexus | II | Dechlorosoma / Dechloromonas | 0.52 |
| CAN_3_bin.203 | JAEULV01 | II | Magnetospirillum | 0.94 |
| CAN_1_bin.77 | UBA5066 | II | Anaeromyxobacter | 0.39\* |
| CAN_3_bin.221 | OLB5 | II | Chryseobacterium / Bacteroidota | 0.92 |
| CAN_5_bin.112, CAN_5_bin.147 | JJ008 / SpSt-398 | II | Flavisolibacter (Bacteroidota) | 0.97 |
| coasm_bin.260 | JAAEKA01 | II | Caldilinea (Chloroflexota) | 1.00 |
| coasm_bin.481 | Leptovillus | II | Rhodothermus / Sphaerobacter | 0.98 |

\* low *sister* support only — Clade II clan membership itself is unambiguous (the
Clade II clan = 233-tip bipartition, 0 Clade I anchors, SH 0.84).

**What this changed vs the two prior artifacts**
- vs **HMM / figure** (`nosz_clades.tsv`, was 1 I / 12 II): agrees on 12/13; the one
  change is CAN_5_bin.40 → Clade I (the gene the HMM itself flagged low-confidence,
  margin −24).
- vs **diamond** (`nosz_clades_updated.tsv` / `partA`, was 7 I / 6 II): the tree
  **rejects 5 of the 7** diamond Clade-I calls (CAN_1_bin.98, coasm_bin.185,
  CAN_2_bin.6, CAN_5_bin.70, CAN_3_bin.203 → all Clade II). Diamond was misled by a
  genus→clade lookup that assumed Rhodocyclaceae nosZ are typical/Clade I; the actual
  Dechloromonas / Dechlorosoma / Accumulibacter / Magnetospirillum reference *sequences*
  fall inside the Clade II clan (with Wolinella, Anaeromyxobacter, Bacteroidota).

**NosZREF-1577D re-examined (Part D step 2):** its core Rhodocyclaceae —
Dechloromonas (3/3), Dechlorosoma (4/4), Accumulibacter (2/2), Magnetospirillum (3/3)
— are all Clade II in the tree, so the `1577D = Clade II` grouping is **correct**, not a
mislabel. (Acidovorax, Alicycliphilus place in Clade I, but they are not 1577D's basis.)

**All three artifacts now carry this single call set** (`clade`/`final_clade`), with each
method's own call preserved in side columns (`clade_hmm`, `clade_diamond`, `tree_clade`).

## Part B — Clade III L-NosZ in all 276 bins
- Searched 276 dereplicated MAGs (935,110 protein ORFs via pyrodigal meta mode)
- He 269NosZ HMM (E ≤ 1×10⁻¹⁰): **180 hits** in 276 bins
- CuZ motif: 126 = DXHH (C-NosZ), 1 = GXHH candidate, 53 = no clear CuZ (fragments)
- **GXHH candidate** (coasm_bin.312::contig_62613_pilon_70, 201 aa):
  - No CuA motif detected → truncated ORF
  - Diamond vs Chee+Orellana: **73.7% identity** to Runella (Bacteroidota clade II)
  - FAILS the <35% identity threshold → **false positive** (clade II C-NosZ fragment)

**Clade III L-NosZ: NOT DETECTED** in any of the 276 dereplicated bins.

## Figure integration
- `gene_rpkm_per_sample_claded.tsv` generated with rows:
  `nirB/D`, `nirS`, `nirK`, `norC` (= max norB/C), `nosZ_cladeI`, `nosZ_cladeII`
- figure `04_plot_heatmap.py` regenerated with per-sample × per-MAG layout
- nosz_clades.tsv unchanged (prior calls confirmed)

## Methods reference
- He et al. 2025, *Nature* 646:152 (L-NosZ / clade III discovery + GraftM refpkg)
- Chee & Orellana 2014 (336 C-NosZ clade I/II reference sequences)
- CuA motif: `C..FC...H.EM`; CuZ: `D.HH` (C-NosZ) vs `G.HH` (L-NosZ)
