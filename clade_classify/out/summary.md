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

## Part A — K00376 nosZ clade I / II reclassification
All 13 K00376 genes were verified against the He 269NosZ HMM and Chee+Orellana
reference set (diamond blastp).

**CuZ motif (definitive N₂OR marker):**
- All 13 = CuA present + CuZ = DXHH → confirmed C-NosZ (not L-NosZ / clade III)

**Clade I / II assignments (from nosz_clades.tsv, validated):**
| MAG | Genus | Clade | Confidence |
|-----|-------|-------|------------|
| coasm_bin.55 | Competibacter | I | low (margin +7.9) |
| CAN_5_bin.40 | Giesbergeria | II | low (margin −24) |
| All others | — | II | high (margin −195 to −366) |

The prior HMM-based calls (nosz_clades.tsv) are correct; diamond-based I/II assignment
was unreliable due to phylogenetic complexity of Betaproteobacteria (Accumulibacter,
Azonexus, Dechloromonas all lie at the clade I/II boundary).

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
