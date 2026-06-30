# Species-novel MAGs of the CAN system — summary

Seven dereplicated MAGs were classified by GTDB-Tk to a named or placeholder genus but carry **no species assignment** (`s__`), with placement *fully defined by tree topology* and no ANI hit to any GTDB reference genome — i.e. each represents a candidate novel species.

### Table N1. Taxonomy, novelty and genome quality

| iterativeID | bin | phylum | genus | compl % | contam % | size (Mb) | contigs | N50 (kb) | GC % | RED |
|---|---|---|---|---|---|---|---|---|---|---|
| Ca_Competibacter.51 | coasm_bin.55 | Pseudomonadota | Competibacter | 99.11 | 1.89 | 4.17 | 3 | 3379 | 53.0 | 0.989 |
| Desulfobacillus.1_m | CAN_2_bin.6 | Pseudomonadota | Desulfobacillus | 99.70 | 0.37 | 3.16 | 4 | 950 | 65.0 | 0.968 |
| Dokdonella.14 | CAN_1_bin.210 | Pseudomonadota | Dokdonella | 99.55 | 0.45 | 3.94 | 13 | 729 | 69.0 | 0.990 |
| Thermoanaerobaculia.1_m | CAN_1_bin.77 | Acidobacteriota | UBA5066 | 97.68 | 0.28 | 5.99 | 4 | 1686 | 70.0 | 0.980 |
| Terrimonas.22 | CAN_5_bin.112 | Bacteroidota | JJ008 | 100.00 | 0.05 | 4.03 | 1 | 4029 | 42.0 | 0.996 |
| Ca_Sarcinithrix.1 | coasm_bin.260 | Chloroflexota | JAAEKA01 | 98.28 | 1.54 | 5.81 | 7 | 3225 | 64.0 | 0.967 |
| Ca_Leptovillus.9 | coasm_bin.481 | Chloroflexota | Leptovillus | 100.00 | 4.49 | 6.26 | 6 | 3282 | 57.0 | 0.993 |

### Table N2. Abundance trajectory and N-cycle inventory

| iterativeID | CAN_1 % | CAN_2 % | CAN_3 % | CAN_4 % | CAN_5 % | trend | N-role | denitrification / DNRA | nosZ |
|---|---|---|---|---|---|---|---|---|---|
| Ca_Competibacter.51 | 0.385 | 0.363 | 0.373 | 1.225 | 3.750 | ↑ | Denitrifier/PAO | nar, nirS, nirB/D(DNRA), nor, nosZ | Clade I |
| Desulfobacillus.1_m | 1.717 | 0.639 | 0.724 | 0.750 | 0.032 | ↓ | Denitrifier/PAO | nap, nirS×2, nirB/D(DNRA), nosZ | Clade II (II.D) |
| Dokdonella.14 | 0.977 | 1.339 | 1.323 | 1.101 | 0.176 | ↓ | Denitrifier/PAO | nirK, nor | — |
| Thermoanaerobaculia.1_m | 1.053 | 1.288 | 1.253 | 1.997 | 4.952 | ↑ | Denitrifier/PAO | nor, nosZ | Clade II (II.G) |
| Terrimonas.22 | 0.030 | 0.011 | 0.015 | 0.759 | 2.777 | ↑ | Denitrifier | nosZ | Clade II (II.C) |
| Ca_Sarcinithrix.1 | 0.008 | 0.009 | 0.010 | 0.867 | 2.178 | ↑ | Denitrifier/PAO | nosZ | Clade II (II.D) |
| Ca_Leptovillus.9 | 2.405 | 1.134 | 1.280 | 3.931 | 4.499 | ↑ | Denitrifier/PAO | nosZ | Clade II (II.D) |

### Table N3. Functional-marker completeness (genes present / genes screened)

| iterativeID | Acetate uptake | Propionate uptake | Glycogen synthesis | Glycogen degradation | PHA cycling | Gluconeogenesis | Polyphosphate | Phosphate transport | Type IV pili | B12 transport | B12 biosynthesis | Biotin biosynthesis | Glycoside hydrolases |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Ca_Competibacter.51 | 3/3 | 3/5 | 3/3 | 3/3 | 2/4 | 3/6 | 1/3 | 6/6 | 11/11 | 4/7 | 4/5 | 4/4 | 1/6 |
| Desulfobacillus.1_m | 3/3 | 2/5 | 3/3 | 3/3 | 3/4 | 3/6 | 1/3 | 6/6 | 11/11 | 3/7 | 3/5 | 4/4 | 0/6 |
| Dokdonella.14 | 3/3 | 2/5 | 0/3 | 0/3 | 4/4 | 4/6 | 2/3 | 5/6 | 11/11 | 4/7 | 0/5 | 3/4 | 1/6 |
| Thermoanaerobaculia.1_m | 2/3 | 2/5 | 2/3 | 3/3 | 1/4 | 5/6 | 2/3 | 5/6 | 8/11 | 4/7 | 0/5 | 2/4 | 1/6 |
| Terrimonas.22 | 2/3 | 1/5 | 3/3 | 2/3 | 1/4 | 4/6 | 1/3 | 0/6 | 0/11 | 3/7 | 0/5 | 0/4 | 3/6 |
| Ca_Sarcinithrix.1 | 1/3 | 2/5 | 2/3 | 3/3 | 1/4 | 3/6 | 2/3 | 5/6 | 0/11 | 0/7 | 0/5 | 1/4 | 1/6 |
| Ca_Leptovillus.9 | 2/3 | 2/5 | 3/3 | 3/3 | 1/4 | 5/6 | 2/3 | 5/6 | 1/11 | 0/7 | 0/5 | 1/4 | 2/6 |

_Functional markers: KOfam (pyhmmer, adaptive thresholds) over Bakta proteomes. Genomic potential only — no metatranscriptome. CAZy/MEROPS/eggNOG deep annotation was not run, so secreted-protease and EPS inventories are out of scope; glycoside-hydrolase counts are KO-level only._
