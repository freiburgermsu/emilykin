# EmilyKin EBPR 16S (V4-V5) ASVs vs the BV-BRC embedding space — findings

## 1. Executive summary

An enhanced-biological-phosphorus-removal (EBPR) activated-sludge bioreactor — a PAO/GAO community dominated by *Candidatus* Accumulibacter, *Ca_*Competibacter, *Thauera*, and *Dechloromonas*-adjacent Rhodocyclaceae — was sequenced for 16S V4-V5 amplicons (515F/926R), yielding **3,950 ASVs** across **77 samples**, each embedded with Nucleotide-Transformer-v2-500m and matched by cosine similarity against a **region-matched** BV-BRC reference space (97,624 unique V4-V5 inserts excised by in-silico PCR). **Sequence placement is highly confident community-wide**: mean best cosine **0.99489**, median **0.99850**, **98.71%** of ASVs match a reference at cosine ≥0.97, **791 ASVs reach exact cosine = 1.0**, and only **1.11% (44 ASVs)** are genuinely novel (best <0.90). **Yet taxonomic-name transfer is weak**: abundance-weighted genus concordance is only **0.395** and family **collapses to 0.156**, because the abundant members are uncultured candidate PAOs/GAOs whose near-exact BV-BRC matches are unclassified MAGs with blank family labels (e.g., the single most abundant ASV — *Ca_*Competibacter, 42.9% — hits an "uncultured Gammaproteobacteria bacterium" MAG at cosine 1.0000). **This INVERTS the prior cultured-methanogen anaerobic-digester study** (codiffusion), where the dominant taxa were cultured *Methanobacterium* and abundance-weighted genus concordance was HIGH (~0.885): in both communities sequence placement is confident (median best cosine ~0.998), but taxonomic-NAME transferability depends on whether the dominant members happen to be cultured and labeled in BV-BRC. The fair signals here are the cosine itself and any-top20 concordance (genus 0.664), both of which confirm the right organisms are present in the neighborhood — only the #1-ranked reference name is missing.

## 2. Methods recap

Each of the 3,950 md5-named ASVs (V4-V5, 515F/926R, ~374 bp expected insert) was embedded with **Nucleotide-Transformer-v2-500m** and scored by cosine similarity against a **region-matched** BV-BRC reference space — 97,624 unique V4-V5 inserts excised from BV-BRC genomes by in-silico PCR with the same 515F/926R primers (the identical embedding mapping reused from the prFBA/codiffusion pipeline). ASVs were also scored against the **full-length (naive)** BV-BRC embeddings for comparison. The study's **MiDAS taxonomy** (`taxonomy.csv`) is the concordance ground truth; concordance comparison **normalizes *Candidatus*/*Ca_* prefixes and NCBI phylum synonyms** before matching MiDAS ranks to best-hit ranks at genus/family/order/class/phylum. Concordance is reported three ways: best-hit (rank-1 reference), any-top20 (concordant if any of the top-20 neighbors matches), and **abundance-weighted**, where abundance = each ASV's `rel_ab` (max % across the 77-sample feature table). The `rel_ab` percent column merged 1:1 onto the concordance table with 0 missing rows.

## 3. Findings

### 3a. Match quality and region-matching

**Match quality (best cosine vs region-matched space, N=3,950):**
- Mean **0.99489**, median **0.99850**, p5 **0.98830** (min 0.7754, max 1.0).
- **%≥0.97 = 98.71%** (3,899 ASVs); **%≥0.99 = 93.01%** (3,674 ASVs).
- Novel fraction (best <0.90) = **1.11%** = **44 ASVs**.

**Region-matching benefit (region-matched vs naive full-length):**
- Naive mean **0.98954** (median 0.9907) → region mean **0.99489**: **mean gain = +0.00535** (median gain +0.0063).
- **94.99%** of ASVs improved (3,752 / 3,950); region ≥ naive for 95.75%; only 4.25% regress.
- **Naive ceiling = 0.998** with **0 ASVs at cosine 1.0**; region-matching breaks that ceiling with **791 ASVs at exact cosine = 1.0**, recovering near-exact V4-V5 inserts that full-length embeddings cannot resolve.

**5 most novel ASVs (lowest best_cosine)** — all MiDAS genus "bacteria", all mapping to the same uncultured hit, and the only systematic exception to the region gain (naive ~0.93–0.94 is *higher* than region ~0.78 because the region-matched space lacks a close Atopobiaceae V4-V5 insert):

| ASV (md5) | best_cosine | naive | best_hit_organism |
|---|---|---|---|
| 9b13e5021d68ef93844f65f09f57923c | 0.7754 | 0.9419 | Atopobiaceae bacterium AF91-08IFCA |
| cd5a90e48ecab9609ecaaeee72c17cfb | 0.7822 | 0.9355 | Atopobiaceae bacterium AF91-08IFCA |
| 00db0c0419e71a3222dc483990c535b9 | 0.7842 | 0.9365 | Atopobiaceae bacterium AF91-08IFCA |
| eef646ce737d72b0aa1a8c0ec4525af2 | 0.7852 | 0.9351 | Atopobiaceae bacterium AF91-08IFCA |
| 76da683c7b87cc132359e85e4d6b665c | 0.7905 | 0.9297 | Atopobiaceae bacterium AF91-08IFCA |

### 3b. Concordance: confident sequences, weak name transfer

**Concordance by rank (`concordance_by_rank.json`, n_asv = 3,950):**

| Rank | n_evaluable | best-hit | any-top20 | abundance-weighted (best) |
|---|---|---|---|---|
| genus | 2,090 | 0.332 | 0.664 | **0.395** |
| family | 3,107 | 0.291 | 0.631 | **0.156** |
| order | 3,566 | 0.304 | 0.538 | 0.162 |
| class | 3,704 | 0.399 | 0.656 | 0.458 |
| phylum | 3,775 | 0.694 | 0.842 | 0.817 |

Placeholder fractions: **genus 0.471**, **family 0.213**.

**Why names fail despite near-perfect cosines** — the drop is concentrated in the most abundant ASVs:
- **1,724** ASVs have best_cosine ≥ 0.999; of these **693 (40.2%) have a blank besthit_family**, and those 693 carry **39.3% of total community abundance** — confident sequence hit, but the matched MAG has no family name to transfer.
- Family abundance-weighted (0.156) is *lower* than its own unweighted best-hit rate (0.291): the abundant taxa actively drag family concordance down.
- The single most abundant ASVs are exactly these failures: *Ca_*Competibacter (rel_ab 42.9%, cos 1.000) → blank-family MAG; the next *Ca_*Competibacter (10.9%, cos 0.997) and *Hydrogenophaga* (9.3%, cos 0.9995) likewise hit blank-family MAGs. Even where the candidate name matches (*Ca_*Accumulibacter → "Candidatus Accumulibacter", genus-concordant), the reference lacks a family — driving the genus 0.395 → family 0.156 cliff.

**Top-20 rescue confirms placement is right, naming is the bottleneck.** Among the 2,090 genus-evaluable ASVs, **694 (33.2%)** are best-discordant but any-top20-concordant — i.e., **49.7%** of all best-discordant evaluable ASVs are rescued by scanning the top-20 neighbors. The correct genus is in the local embedding neighborhood; it simply isn't the #1 reference (cosine ties among near-identical MAGs, or the #1 reference is an unnamed MAG). Phylum-level abundance-weighted concordance stays high (**0.817**) because phylum names survive even for uncultured MAGs, isolating the failure to genus/family naming rather than sequence biology.

### 3c. The dominant PAO/GAO guild

Top 15 ASVs by max relative abundance (% across 77 samples). Region-matched and naive cosines are **identical for all 15**, so the V4-V5 excision did not degrade placement of the dominant members. `genus_ok`/`family_ok` use the any20 criterion.

| rank | rel_ab% | midas_genus | best_hit | cosine | genus_ok(any20) | family_ok(any20) |
|---|---|---|---|---|---|---|
| 1 | 42.91 | Ca_Competibacter (GAO) | uncultured Gammaproteobacteria bacterium M… | 1.0000 | yes | yes |
| 2 | 33.27 | Thauera | Thauera mechernichensis #50 (MBG-DUTH)/TV-… | 1.0000 | yes | no |
| 3 | 21.05 | midas_g_31688 | Roseibium sp. CAU 1637 | 0.9985 | no (n/a) | no |
| 4 | 16.67 | Ca_Accumulibacter (PAO) | Candidatus Accumulibacter meliphilus UW14 | 1.0000 | yes | yes |
| 5 | 12.36 | Pseudoxanthomonas | Pseudoxanthomonas mexicana isolate_71 | 1.0000 | yes | no |
| 6 | 10.90 | Ca_Competibacter (GAO) | uncultured bacterium CTOTU24409 | 0.9971 | yes | yes |
| 7 | 9.77 | Ca_Accumulibacter (PAO) | Candidatus Accumulibacter sp. UW20 | 1.0000 | yes | yes |
| 8 | 9.30 | Hydrogenophaga | Betaproteobacteria bacterium FK_Sedi_B_Bin… | 0.9995 | yes | yes |
| 9 | 8.97 | midas_g_171 | Ignavibacteria bacterium AS-MICRO-OLR.79 | 1.0000 | no (n/a) | no (n/a) |
| 10 | 8.39 | Aeromonas | Aeromonas hydrophila strain 3019 | 1.0000 | yes | yes |
| 11 | 8.09 | Ca_Competibacter (GAO) | uncultured Gammaproteobacteria bacterium M… | 0.9966 | yes | yes |
| 12 | 6.93 | midas_g_71310 | Patescibacteria group bacterium UBA5532 | 0.9990 | no (n/a) | no (n/a) |
| 13 | 6.69 | Ca_Leptovillus | Chloroflexi bacterium strain Kalu_18-Q3-R1… | 1.0000 | no | no |
| 14 | 6.04 | Zoogloea | Proteobacteria bacterium strain SZAS-76 | 0.9995 | yes | yes |
| 15 | 5.62 | Azonexus | uncultured bacterium CTOTU26046 | 0.9990 | yes | no |

**Functional-guild readout:**
- **Every guild member is recognized by SEQUENCE.** All named EBPR-guild taxa sit at cosine ≈ 1.000 (*Ca_*Competibacter 1.0000/0.9971/0.9966; *Ca_*Accumulibacter 1.0000/1.0000; *Thauera* 1.0000; Rhodocyclaceae relatives *Zoogloea*/*Azonexus* 0.9990–0.9995).
- ***Ca_*Accumulibacter (PAO) is the standout success:** both ASVs (ranks 4, 7) hit named *Candidatus Accumulibacter* genomes (UW14 meliphilus, UW20) — genus NAME transfers on the best hit too; family fails only because the matching Accumulibacter MAGs carry no family label (rescued by another top-20 hit).
- ***Ca_*Competibacter (GAO), the single most abundant taxon (42.9%), is recognized but UN-NAMED:** all three Competibacter ASVs (ranks 1, 6, 11) match "uncultured Gammaproteobacteria/bacterium" MAGs with empty besthit_genus/family; concordance is rescued only via any20. This is the core PAO/GAO contrast — confident sequence, no name on the nearest genome.
- ***Thauera* (rank 2):** genus name transfers, but family does not — BV-BRC files it under Zoogloeaceae vs MiDAS Rhodocyclaceae, a genuine reference-taxonomy disagreement, not an embedding failure.
- **Cultured isolates behave like the digester study:** *Aeromonas hydrophila* (rank 10), *Pseudoxanthomonas mexicana* (5), *Thauera* (2) are exactly where the name transfers at genus level, mirroring the codiffusion methanogen result.
- **3 of 15 (ranks 3, 9, 12) are MiDAS placeholder genera** (midas_g_31688/171/71310) — not genus-evaluable at all despite cosine ≥ 0.9985, contributing to the low abundance-weighted genus concordance.

Bottom line: of the top-15 abundant ASVs, **13/15 are placed at cosine ≥ 0.997** and **10/15 reach genus-level NAME concordance via any20**, but the dominant GAO (*Ca_*Competibacter, 42.9%) and several abundant candidate taxa map to label-less MAGs, so the name fails to transfer exactly where abundance is highest.

### 3d. BV-BRC coverage gaps for uncultured EBPR taxa

**Low-cosine "true gap" ASVs.** Only **44 / 3,950 have best_cosine < 0.95, and all 44 are also < 0.90** — there are **zero ASVs in the 0.90–0.95 band**, a bimodal distribution (tight high-similarity mass + a small far-off cluster). The 8 lowest-cosine ASVs are one homogeneous group — all length **543 bp**, all MiDAS-classified only to "bacteria", all matching the same *Atopobiaceae bacterium AF91-08IFCA* MAG. The 543 bp length (vs ~374 bp expected V4-V5 insert) flags these as off-target/chimeric or non-target amplicons; the embedding correctly fails to find a close region-matched reference, so these are genuine reference-space voids, not mis-mappings.

**High-confidence placement with no taxonomic name (the core gap).** Of the **3,107 family-evaluable ASVs, 3,012 have best_cosine ≥ 0.99**. Of those, **1,095 (36.4%) have NO family-concordant hit anywhere in the top-20** — confident placement, no transferable family name. **688 of those 1,095 (62.8%)** map to an uncultured/MAG best hit (679 organism names contain "uncultured"/"bacterium"/"archaeon"/"metagenome"; 514 have a blank besthit_family). This is the BV-BRC coverage gap made quantitative: high sequence confidence collides with poorly-classified candidate-PAO/GAO MAGs.

**Low-cosine ASVs are completely enriched for placeholder genera.** With placeholder = MiDAS `midas_g_*`/blank/unranked "bacteria"/"archaea" (overall placeholder-genus fraction = **0.358**):
- best_cosine < 0.95: placeholder fraction = **1.000** (44/44)
- best_cosine < 0.90: placeholder fraction = **1.000** (44/44)
- best_cosine ≥ 0.95: 0.351; ≥ 0.99: 0.329

Enrichment = **2.79× over background**; Fisher exact OR = ∞ (44 placeholder / 0 named below 0.95), **p = 1.5e-20**. The rare true sequence-space gaps and the placeholder-taxonomy gaps coincide perfectly: when no good reference exists, MiDAS also failed to name it. Net picture: only ~1.1% of ASVs are genuine reference voids, all MiDAS-unnamed off-length amplicons — the dominant coverage gap is **name absence, not sequence absence**.

## 4. Interpretation & caveats

**Sequence placement is confident for both communities.** Whether the dominant members are cultured methanogens (digester) or uncultured candidate PAOs/GAOs (this EBPR study), the NT-v2 embedding places V4-V5 ASVs on a near-exact region-matched reference (median best cosine ~0.998; here 98.7% ≥ 0.97, 791 exact-1.0 matches). The biology — "which sequence is this" — is solved robustly in both cases.

**Taxonomic-NAME transfer is weak here, and it INVERTS relative to the digester.** In the codiffusion digester, the abundant taxa were cultured *Methanobacterium*, so high-cosine matches carried transferable cultured names and abundance-weighted genus concordance was HIGH (~0.885). In this EBPR community the abundant taxa (*Ca_*Competibacter 42.9%, *Ca_*Accumulibacter, *Thauera*, etc.) are largely uncultured candidate organisms whose near-exact BV-BRC matches are **unclassified MAGs with blank family labels** (≈54% abundance-weighted blank-family among cos≥0.99 ASVs) or carry *Candidatus*/name-format mismatches. The result is the inverse: abundance-weighted genus concordance **0.395**, family **0.156** — the abundant members, despite cosine ~1.0, drag concordance down rather than up.

**Report the fair signals, not the best-hit name alone.** Best-hit genus/family concordance understates the embedding's success because the #1 reference is often an unnamed MAG. The honest readouts are (i) **the cosine itself** (median 0.9985 — the sequence match is real) and (ii) **any-top20 concordance** (genus 0.664, family 0.631 — the correct name is usually present in the neighborhood, rescuing 49.7% of best-discordant genus-evaluable ASVs). Phylum-level transfer remains high (abundance-weighted 0.817) because phylum names survive for uncultured MAGs.

**Caveats.** (i) The 44 genuine sequence-void ASVs are off-length (543 bp) non-target/chimeric amplicons, not informative biology. (ii) MiDAS placeholder genera (overall 0.471 of best hits at genus, 0.358 of ASVs by MiDAS label) mean a large block of ASVs is **not genus-evaluable on either side**, so concordance rates are computed only over evaluable subsets. (iii) Some family discordances (e.g., *Thauera* Zoogloeaceae vs Rhodocyclaceae) are reference-taxonomy disagreements between BV-BRC and MiDAS, not embedding errors. (iv) Abundance is the per-ASV max % across the 77 samples, so abundance-weighting reflects peak community dominance.

## 5. Outputs

- `/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/asv_top20_hits.json` — per-ASV top-20 region-matched neighbors (cosine, organism, genome, taxon_id, n_genomes, ref_seq_len)
- `/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/asv_summary.csv` — per-ASV best hit, best/naive cosine, genus_concordant, novel flag
- `/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/asv_concordance.csv` — per-ASV abundance + MiDAS/best-hit ranks + best/any20/evaluable flags per rank (Ca_-normalized)
- `/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/concordance_by_rank.json` — best/any20/abundance-weighted concordance + placeholder fractions per rank
- `/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/findings_stats.json` — match-quality stats (best_cosine percentiles, naive-vs-region gain, match rates)

---

# Appendix: adversarial verification

All recomputations are complete and consistent. Here is the verification.

## Verification

**Data:** 3,950 ASVs, V4-V5 (515F/926R), 97,624 region-matched BV-BRC refs. `rel_ab` merged from `taxonomy.csv` onto `asv_concordance.csv` with **0 missing** (clean 1:1 join).

**Claim 1 — median best_cosine ≥0.998 and ≥98% match at cosine≥0.97 — CONFIRMED.**
Recomputed median best_cosine = **0.9985** (≥0.998 ✓). Fraction at cosine ≥0.97 = **0.9871** (3899/3950, ≥98% ✓). (For reference, ≥0.99 = 0.9301.)

**Claim 2 — abundance-weighted genus ~0.40, family ~0.16, both LOW vs single-best-hit — CONFIRMED.**
Weighting `<rank>_best_concordant` by `rel_ab` over evaluable ASVs:
- genus: AW-best = **0.3946** (n_eval=2090); unweighted single-best-hit = 0.3316.
- family: AW-best = **0.1563** (n_eval=3107); unweighted = 0.2910.
Family AW is *lower* than its own unweighted rate — the abundant taxa drag family concordance down, confirming the contrast. (order AW=0.1615, class=0.4581, phylum=0.8167.) Matches `concordance_by_rank.json` exactly.

**Claim 3 — dominant Ca_Accumulibacter ASV(s) match Accumulibacter ref at cosine ~1.0 — CONFIRMED.**
The two dominant Accumulibacter ASVs are community abundance ranks **#4 and #7**:
- `8b541735…` rel_ab=16.67% → *Candidatus Accumulibacter meliphilus UW14*, cosine = **1.0000**.
- `91a3df65…` rel_ab=9.77% → *Candidatus Accumulibacter sp. UW20*, cosine = **1.0000**.
All top-5 Accumulibacter ASVs hit *Candidatus Accumulibacter* references at cosine ≥0.9995 (one minor exception maps to Propionivibrio, a known Accumulibacter relative).

**Claim 4 — high-cosine (≥0.99) abundant ASVs map to uncultured/MAG refs with blank besthit_family — CONFIRMED.**
Of 3,674 ASVs at cosine ≥0.99, **46.2%** (1699) have a blank `besthit_family`; **abundance-weighted = 54.4%**. Among the 50 most-abundant cos≥0.99 ASVs, **27/50 (54%)** have blank family. Blank `besthit_genus` = 59.6% (AW 56.7%). The "large share" is verified, strongest when abundance-weighted. (E.g., the #1 overall ASV, Ca_Competibacter rel_ab=42.9%, hits "uncultured Gammaproteobacteria bacterium" at cosine 1.0000.)

**Claim 5 — region-matching improves cosine over naive for the large majority — CONFIRMED.**
Region > naive for **94.99%** of ASVs (3752/3950); region ≥ naive for 95.75%; only 4.25% regress. Mean gain = **+0.0053**, median gain = +0.0063. Matches `findings_stats.json` (region_match_cosine_gain_mean = 0.0053).

**Overall:** 5/5 CONFIRMED. The core narrative holds: near-exact sequence placement (median cosine 0.9985, 98.7% ≥0.97) coexists with low abundance-weighted taxonomic-name transfer (genus 0.39, family 0.16) because the dominant uncultured PAO/GAO ASVs — though matched at cosine ~1.0 (e.g., Accumulibacter, Competibacter) — land on poorly-classified BV-BRC MAGs (≈54% abundance-weighted blank-family among cos≥0.99).
