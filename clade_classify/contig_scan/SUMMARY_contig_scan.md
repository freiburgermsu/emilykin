# Part C — clade III L-NosZ search over ALL assembly contigs (incl. unbinned)

Higher-resolution extension of Part B. Part B searched only the **276 dereplicated
MAGs** (binned + dereplicated fraction). This search covers the **entire assembly**
— the co-assembly plus all five per-sample assemblies — so it also sees the
**unbinned contigs** and any bins removed during dereplication, which is exactly
where a novel, low-abundance, poorly-binned L-NosZ host would hide.

Run: `~/Documents/py_venv/bin/python clade_classify/scan_all_contigs.py` (≈ 8.5 min,
48 threads). Outputs: `partC_motif.tsv`, `partC_lnosz_hits.tsv`, `scan.log`, and the
gitignored `all_contigs_orf.faa`.

## Method (identical to Part B, for direct comparability)

pyrodigal meta ORFs → He 269NosZ HMM (pyhmmer `hmmsearch`, E ≤ 1×10⁻¹⁰) →
CuA/CuZ motif (`DXHH` = C-NosZ clade I/II ; `GXHH` = L-NosZ clade III) →
each `GXHH` candidate confirmed by % identity to Chee+Orellana C-NosZ (< 35 % = L-NosZ).
The diamond confirm step of Part B is replaced by a Biopython local alignment
(BLOSUM62) — diamond is not installed locally and the < 35 % vs 60–90 % gap is robust
to the exact identity definition.

## Inputs

| Source | File | Contigs |
|--------|------|--------:|
| co-assembly | `meta/results/08a_map_shortreads_co/all_co_contigs.fasta` | 49,089 |
| CAN_1 | `meta/results/08_map_shortreads/CAN_1/CAN_1_contigs.fasta` | 22,932 |
| CAN_2 | `meta/results/08_map_shortreads/CAN_2/CAN_2_contigs.fasta` | 20,952 |
| CAN_3 | `meta/results/08_map_shortreads/CAN_3/CAN_3_contigs.fasta` | 25,254 |
| CAN_4 | `meta/results/08_map_shortreads/CAN_4/CAN_4_contigs.fasta` | 21,570 |
| CAN_5 | `meta/results/08_map_shortreads/CAN_5/CAN_5_contigs.fasta` | 16,510 |
| **total** | | **156,307** |

## Results

| Metric | Part B (276 MAGs) | **Part C (all contigs)** |
|--------|------------------:|-------------------------:|
| ORFs searched | 935,110 | **4,808,027** |
| He-HMM hits (E ≤ 1e-10) | 180 | **909** |
| `DXHH` (C-NosZ, clade I/II) | 126 | **614** |
| `GXHH` (L-NosZ candidate) | 1 | **1** |
| no clear CuZ (fragments) | 53 | **294** (217 no-CuA + 77 CuA-only) |
| **confirmed clade III L-NosZ (< 35 % id)** | **0** | **0** |

≈ 5.1× more ORFs and 5× more N₂OR hits than the MAG-only search — i.e. the unbinned
fraction was genuinely added — yet still **zero** clade III L-NosZ.

### The lone GXHH candidate is again a *Runella* C-NosZ false positive

`CAN_3::ctg17630_19` (219 aa, **no CuA** → truncated ORF) best-matches
*Runella slithyformis* (Bacteroidota, clade II C-NosZ) at **71.3 % identity** —
far above the 35 % L-NosZ ceiling → **C-NosZ fragment, not L-NosZ**. This mirrors
Part B's single false positive (`coasm_bin.312::contig_62613_pilon_70`, 73.7 % to
*Runella*): across two independent searches the *only* thing that ever trips the
`GXHH` regex is a *Runella*-type clade II C-NosZ fragment, not a real L-NosZ.

## Conclusion

**Clade III L-NosZ is NOT present in this system**, now confirmed at full-assembly
resolution (incl. unbinned contigs), not just within the dereplicated MAGs. N₂O
reduction here is exclusively C-NosZ (clade I/II, K00376).

## Remaining (unsearched) sensitivity tiers

This was an ORF-level (pyrodigal meta) HMM search — the same as Part B. The pipeline
*does* call and detect partial ORFs (294 fragment hits prove the HMM is firing on
fragments), so a real L-NosZ would almost certainly have produced at least a `GXHH`
hit. If even higher sensitivity is wanted:
1. **6-frame translated search** (`hmmsearch --dnax` or manual 6-frame → pyhmmer) of
   the same contigs — marginally more sensitive on heavily fragmented contigs.
2. **Raw long reads** (`meta/longreads/*.fastq.gz`, ~118 GB) — the only fraction never
   searched; last resort (½–3 days, translated search of ~150–200 Gbp).

Given the consistent *Runella*-only `GXHH` pattern across the bins and the full
assembly, a missed L-NosZ is unlikely; tiers 1–2 are optional confirmation.
