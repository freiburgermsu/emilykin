# nosZ Clade I / II — unresolved disagreement between the result artifacts

**Status: not yet resolved.** Three shipped artifacts give two different answers
to "how many of the 13 K00376 nosZ genes are Clade I?" This document states the
conflict explicitly, because nowhere else does. (`out/summary.md` lines 29–31
assert a *conclusion* — "HMM calls are correct, diamond was unreliable" — but never
shows that the two `out/` spreadsheets still carry the opposite calls, which 6
genes differ, or why.)

## TL;DR

| Artifact | Method | Clade I | Clade II |
|---|---|---:|---:|
| `gene_ab_figure/data/nosz_clades.tsv` → **what the figure plots** | nucleotide HMM bit-score (`nhmmer`, best clade-I vs best clade-II subclade) | **1** | 12 |
| `clade_classify/out/partA_clade_I_II.tsv` (`final_clade`) | diamond blastp top hit vs Chee+Orellana, taxonomy-informed | **7** | 6 |
| `clade_classify/out/nosz_clades_updated.tsv` (`clade`) | same as partA | **7** | 6 |

So **the figure says 1 Clade I; the two spreadsheets in `out/` say 7.** They
disagree on **6 genes**, all flipping HMM-Clade-II → diamond-Clade-I.

## Per-gene breakdown (all 13 K00376 genes)

HMM = figure call (`nosz_clades.tsv`); diamond = `partA_clade_I_II.tsv` `final_clade`.

| Gene | Genus | HMM (figure) | HMM margin | HMM subclade | Diamond top hit (Chee+Orellana) | %id | Diamond call | Agree? |
|---|---|:---:|---:|:---:|---|---:|:---:|:---:|
| coasm_bin.55 | *Competibacter* | **I** | +7.9 | II.G | *Rhodoferax ferrireducens* | 74 | **I** | ✅ both I |
| CAN_1_bin.77 | UBA5066 | II | −311.6 | II.G | *Anaeromyxobacter* sp. | 81 | II | ✅ both II |
| CAN_3_bin.221 | OLB5 | II | −246.6 | II.B | *Hydrogenobacter thermophilus* | 69 | II | ✅ both II |
| CAN_5_bin.112 | JJ008 | II | −266.3 | II.C | *Dyadobacter fermentans* | 74 | (unknown)→II | ✅ both II |
| CAN_5_bin.147 | SpSt-398 | II | −276.9 | II.C | *Haliscomenobacter hydrossis* | 74 | (unknown)→II | ✅ both II |
| coasm_bin.260 | JAAEKA01 | II | −195.5 | II.D | *Caldilinea aerophila* | 76 | II | ✅ both II |
| coasm_bin.481 | *Leptovillus* | II | −213.7 | II.D | *Rhodothermus marinus* | 60 | II | ✅ both II |
| **CAN_1_bin.98** | ***Accumulibacter*** | II | −364.4 | **II.D** | *Ca.* Accumulibacter phosphatis | 92 | **I** | ❌ **FLIP** |
| **coasm_bin.185** | ***Accumulibacter*** | II | −358.2 | **II.D** | *Ca.* Accumulibacter phosphatis | 91 | **I** | ❌ **FLIP** |
| **CAN_2_bin.6** | *Desulfobacillus* | II | −365.6 | **II.D** | *Dechloromonas aromatica* RCB | 86 | **I** | ❌ **FLIP** |
| **CAN_5_bin.70** | *Azonexus* | II | −340.5 | **II.D** | *Dechloromonas aromatica* RCB | 90 | **I** | ❌ **FLIP** |
| **CAN_3_bin.203** | JAEULV01 | II | −274.2 | **II.D** | *Magnetospirillum magneticum* | 70 | **I** | ❌ **FLIP** |
| **CAN_5_bin.40** | *Giesbergeria* | II *(low)* | −24.2 | II.G | *Alicycliphilus denitrificans* | 87 | **I** | ❌ **FLIP** |

7 genes agree (1 as I, 6 as II); **6 genes flip** (II → I).

## Why they disagree — the most likely root cause

**The conflict is concentrated on the reference "II.D" subclade.** Five of the six
flips have HMM best-match `II.D` (the 6th, *Giesbergeria*, is a borderline II.G at
margin −24). The `II.D` reference file (`NosZREF-1577D`, *Dechloromonas /
Dechlorosoma*) was **assumed** to be Clade II in the original gene-abundance step —
not by an explicit key, but by inference from file naming and taxon composition
(documented in `gene_ab_figure/METHODS_nosz_clades.md` §1: *"Clade grouping — an
interpretation, not an edit"*).

But *Dechloromonas*, *Dechlorosoma*, *Accumulibacter*, *Azonexus*, and
*Alicycliphilus* are **Rhodocyclaceae/Comamonadaceae denitrifiers that carry
canonical, "typical" Clade I nosZ** (Sanford et al. 2012; Jones et al. 2013).
Diamond — which compares against *taxonomically labeled* reference proteins at
86–92% identity — recovers exactly that. So the disagreement is plausibly not
"diamond is unreliable" but **"the reference's II.D group is mislabeled as Clade
II when it is really Clade I,"** which would make the HMM (and therefore the
figure) systematically miscall these genes as Clade II.

In other words: the HMM gave huge, confident-looking margins toward II.D
(−274 to −365 bits), but a confident match to a *possibly-mislabeled reference
group* is confidently wrong, not reliable. `summary.md`'s claim that the HMM
calls are "correct" was asserted, not demonstrated — no tree-based adjudication
is shown, even though a query-placed tree was built (`out/nosz_tree.nwk`,
`out/nosz_aligned.afa`).

## Naming-collision trap (flag for collaborators)

The diamond top hit for both *Accumulibacter* genes is
`Candidatus_Accumulibacter_phosphatis_clade_IIA_str._UW-1`. Here **"clade IIA" is
*Accumulibacter*'s own PAO lineage designation, not a nosZ clade.** Do **not**
read it as evidence for nosZ Clade II — *Accumulibacter* nosZ is Clade I.

## Impact

- The **figure** currently shows essentially no Clade I denitrifiers (1 of 13),
  which would read as "this system's N₂O reduction is almost entirely atypical
  Clade II." If the diamond/taxonomy calls are right, the true split is ~7 I / 6
  II, a materially different ecological story (canonical denitrifier-type N₂O
  reduction is well represented, including in the *Accumulibacter* PAOs).
- The **figure, `summary.md`, and `partA_clade_I_II.tsv` contradict each other.**
  Anyone opening the spreadsheet sees 7 Clade I; anyone reading the figure sees 1.

## Recommended resolution (definitive)

Neither best-HMM-score nor best-diamond-hit is authoritative on its own. Settle it
phylogenetically:

1. **Tree-based placement** — place the 13 query proteins in a tree with reference
   sequences whose Clade I / II membership is *independently* established (not the
   `1577*` file-name assumption). Material already exists: `out/nosz_aligned.afa`,
   `out/nosz_tree.nwk`. Assign each query to the Clade I or Clade II monophyletic
   group it falls within. (This is the "placement-preferred" path in
   `PLAN_clade_classification.md` §5 that was built but not used for the final call.)
2. **Re-examine the reference grouping** — verify whether `NosZREF-1577D`
   (*Dechloromonas / Dechlorosoma*) is genuinely Clade II or actually Clade I. If
   the latter, rebuild the Clade I HMM to include these typical-nosZ
   Betaproteobacteria and re-score.
3. **Then reconcile all three artifacts** to one call set and regenerate the
   figure + `summary.md` so they agree.

Until step 3, treat the Clade I count as **unresolved (between 1 and 7)** in any
external communication.
