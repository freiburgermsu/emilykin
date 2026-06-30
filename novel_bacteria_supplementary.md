# Supplementary Information — Novel bacteria of the CAN system

*Functional and ecological characterisation of the seven species-novel MAGs,
applying the analytical framework of the reference supplementary information
(Wang, Gao & Wells, "Integrated Omics Analyses Reveal Differential Gene
Expression and Potential for Cooperation Between Denitrifying Polyphosphate and
Glycogen Accumulating Organisms"; emi15486-sup-0001-supinfo.docx) to this
project's novel organisms.*

Companion data: **Table N1–N3** (`gene_ab_figure/data/novel_bacteria_summary.{tsv,md}`),
**Figure N1** (`gene_ab_figure/novel_denitrification_bubble.pdf`; the metagenomic
analog of reference Figure S3) and **Figure N2**
(`gene_ab_figure/novel_metabolic_potential.pdf`; the analog of reference Figure S4).

> **Scope and caveats.** All statements below describe *metagenomic gene
> potential*, not expression: this study has no metatranscriptome, so the
> reference's transcript-proportion arguments (e.g. "~88 % of *vpr* transcripts
> mapped to CH7") have no analog here. Functional inventories come from KOfam
> (pyhmmer, adaptive per-KO thresholds; HMMER3 v0.12.1) over the Bakta proteomes,
> using a curated marker-KO panel spanning the reference's functional axes; the
> three gluconeogenesis markers reproduce the independent prior KOfam run exactly
> (`kofam_new/new_ko_copynumber.tsv`). Deep CAZy/MEROPS/eggNOG annotation (DRAM)
> was **not** run for this dataset, so secreted-protease and EPS-synthesis
> inventories — central to the reference's CH7/PR6 scavenging argument — cannot be
> resolved here; glycoside-hydrolase rows in Figure N2 are KO-level markers only.

## S1. Seven candidate novel species span four phyla

From the 276 dereplicated MAGs, seven were classified by GTDB-Tk to a named or
placeholder genus but carry **no species assignment** (`s__`), with placement
*fully defined by tree topology* and no average-nucleotide-identity (ANI) hit to
any GTDB reference genome (Table N1; closest-genome ANI = N/A, RED 0.967–0.996).
Each therefore represents a candidate novel species. The set spans four phyla —
three Pseudomonadota (Gammaproteobacteria), one Acidobacteriota, two Chloroflexota
and one Bacteroidota — and all are high-quality genomes (CheckM2 completeness
97.7–100 %, contamination ≤ 4.5 %; *Terrimonas.22* is a single closed 4.03 Mb
contig). Six of the seven are classified Denitrifier/PAO by genome content and one
(*Terrimonas.22*) Denitrifier; none is a flanking heterotroph by the manuscript's
KO logic. Five of the seven increase in relative abundance toward the end of the
time series (CAN_4–CAN_5), echoing the community-level late rise of *nosZ* Clade II
reported in the main text.

## S2. Carbon, energy and storage metabolism

As in the reference, the novel populations are unified by a complete capacity for
short-chain fatty-acid activation and storage-polymer cycling — the metabolic core
of polyphosphate- and glycogen-accumulating organisms. Acetate activation (*ackA–
pta* and/or *acs*) is present in all seven; propionate routing through propionyl-CoA
synthetase (*prpE*) and the methylcitrate/methylmalonyl-CoA branches is present in
six. Glycogen synthesis (*glgC/glgA/glgB*) and degradation (*glgP/glgX/malQ*) and
PHA cycling (*phaA/phaB/phaC*) are widespread (Figure N2). Critically, the
glycogen-accumulating gluconeogenic signature highlighted in reference Table S3 —
NAD-dependent glyceraldehyde-3-phosphate dehydrogenase (*gapA*), PEP carboxykinase
(*pckA*) and pyruvate-phosphate dikinase (*ppdk*) — is recovered across the set:
all seven carry *gapA* + *pckA*, and six also carry *ppdk* (absent only in
*Dokdonella.14*, which instead encodes *pps* and *fbp*). Polyphosphate kinase
(*ppk1* and/or *ppk2*) is present in all seven, while the high-affinity phosphate
transporter *pstSCAB* with the *phoB/phoR* regulon is present in all but
*Terrimonas.22*, supporting the Denitrifier/PAO assignment of the other six.
*Terrimonas.22* lacks the *pst* operon and *ppk1* (retaining only *ppk2*),
consistent with a non-PAO lifestyle.

## S3. Three functional guilds among the novel bacteria

Overlaying the metabolic potential (Figure N2) on the denitrification inventory
(Figure N1) resolves the seven novel MAGs into three ecological guilds that
parallel the organism types of the reference community.

**(i) Biofilm-forming, carbon-cycling denitrifiers (the PAO/GAO-flanking guild).**
The three novel Pseudomonadota (*Ca_Competibacter.51*, *Desulfobacillus.1_m*,
*Dokdonella.14*) and the Acidobacterium (*Thermoanaerobaculia.1_m*) each encode an
essentially **complete type IV pilus biosynthesis operon** (*pilA, pilB, pilC,
pilD, pilM–Q, pilT, pilU*; *Thermoanaerobaculia.1_m* is near-complete, missing
*pilA*, *pilP* and *pilU*). In
the reference, type IV pili (PilQ, PilY1) expressed by the Accumulibacter PAOs and
the flanking population PR6 were proposed to promote granule aggregation and
biofilm formation; the same machinery is genomically present in this guild and,
by the same logic, may contribute to the CAN granule architecture (a hypothesis
that, as in the reference, awaits direct verification). Within this guild
*Ca_Competibacter.51* (coasm_bin.55) is the direct analog of the reference GAO1: a
Competibacteraceae glycogen-accumulating organism carrying the full *gapA/pckA/
ppdk* gluconeogenic route, the most complete denitrification pathway of any novel
MAG (*nar* + *nirS* + *nirB/D* + cNOR + *nosZ*), and — uniquely in this study —
the sole **Clade-I** *nosZ*. It blooms an order of magnitude over the series
(0.38 % → 3.75 %).

**(ii) Filamentous Chloroflexota nitrous-oxide sinks.** *Ca_Leptovillus.9* and
*Ca_Sarcinithrix.1* (Anaerolineae) retain glycogen and gluconeogenic capacity but
encode **only the terminal denitrification step**, a Clade-II.D *nosZ*, and no
upstream nitrate/nitrite/NO-reductase (Figure N1). They notably **lack the type IV
pilus operon** (only the multifunctional prepilin peptidase *pilD* in
*Ca_Leptovillus.9*) and the canonical *btuB/tonB* B12-uptake markers. Both bloom
late (*Ca_Sarcinithrix.1* rises >250-fold, 0.008 % → 2.18 %), positioning the
filamentous Chloroflexota as dedicated late-stage N₂O consumers.

**(iii) A Bacteroidota polysaccharide scavenger.** *Terrimonas.22*
(Chitinophagaceae) is the closest analog of the reference scavengers CH7/PR6: it
encodes **only a Clade-II.C *nosZ*** for N-cycling (a pure N₂O sink), lacks all PAO
phosphate machinery, and is the only novel MAG enriched in secreted-type glycoside
hydrolases — cellulase (*celA*), xylan 1,4-β-xylosidase (*xynB*) and chitinase
(*chiA*) — suggesting a niche in degrading polysaccharides (decaying biomass, EPS,
chitin) rather than competing for primary VFAs. It blooms from 0.03 % to 2.78 %
and is recovered as a single closed contig.

## S4. Putative metabolite exchange (cooperation hypotheses)

The reference frames its flanking populations around metabolite exchange with the
PAOs (B12 provision, EPS/protein scavenging, micronutrient import). Two analogous,
explicitly hypothetical, exchange axes are visible in the novel-bacteria genomes:

- **Vitamin-B12 cross-feeding.** Most novel MAGs encode the TonB-dependent B12
  uptake apparatus (*btuB* and/or *tonB/exbB/exbD*) but **lack the de novo
  cobalamin pathway**; only *Ca_Competibacter.51* (*cobQ/cobS/cobT/cobU*) and
  *Desulfobacillus.1_m* (*cobS/cobT/cobU*) carry biosynthetic *cob* markers. The
  novel community is therefore largely **B12-auxotrophic**, dependent on cobalamin
  produced by a few members (or by the wider community) — the same uptake-versus-
  synthesis asymmetry the reference invokes for its *btuB*-expressing flanking
  populations. (The two Chloroflexota encode neither uptake nor synthesis markers,
  implying either B12-independence or an uncaptured transporter.)
- **Polysaccharide turnover.** The glycoside-hydrolase repertoire of
  *Terrimonas.22* could liberate oligosaccharides from EPS and detrital
  polysaccharide, a scavenging role analogous to the reference's CH7/PR6
  glycoside-hydrolase/protease activity. Confirming this would require the CAZy/
  MEROPS annotation that is unavailable for this dataset.

Cofactor self-sufficiency varies in parallel: the four Proteobacterial/Acidobacterial
guild-(i) members encode multiple biotin-biosynthesis steps (subsets of
*bioB/bioA/bioD/bioF*; all four in *Ca_Competibacter.51* and *Desulfobacillus.1_m*),
whereas the Chloroflexota carry only *bioF* and *Terrimonas.22* none, reinforcing
the picture of streamlined, exchange-dependent novel populations alongside more
self-sufficient carbon-cycling denitrifiers.

## S5. Summary

Applying the reference supplementary framework to the seven species-novel MAGs
recovers the same organising logic — carbon/energy and storage metabolism,
polyphosphate/phosphate handling, type IV pili and biofilm formation, vitamin-B12
exchange, and a graded division between self-sufficient carbon-cyclers and
exchange-dependent scavengers — now mapped onto a denitrification-anchored
community. The novel bacteria are not metabolic bystanders: one is a GAO1-like
glycogen-accumulating denitrifier carrying the study's only Clade-I *nosZ*, and the
filamentous Chloroflexota and the Bacteroidota scavenger together constitute a
late-blooming guild of dedicated nitrous-oxide sinks.
