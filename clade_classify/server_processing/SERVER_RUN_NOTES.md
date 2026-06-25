# Server run notes — nosZ RPKM + depth-copy-number for ALL nosZ bins

Executed on poplar. This records how `compute_nosz_rpkm.py` was actually run and
one correction to the premise in `README_server_processing.md`.

## Premise correction: the focal BAMs did not cover all nosZ bins

`README_server_processing.md` (§0) points the script at the
`aligned_<sample>_sorted.bam` files that `gene_ab_figure/03_count_and_rpkm.py` /
`05_gene_copy_number.py` used. Those BAMs are mapped to the **21-MAG focal
reference** (`gene_ab_figure/data/combined_ref.fa`, 267 contigs), not the full
dereplicated catalogue. The all-bins BED (`allbins_nosz_loci_clade.bed`) spans
**122 MAGs / 136 nosZ loci**; only **13** of those MAGs are in the focal set, so
a first run against the focal BAMs found only **22/136** loci and counted 0 for
the other 114 (the naming check warned about this — it was a *scope* mismatch,
not a contig-naming-style mismatch; both use `<mag_underscores>::<contig>`).

## What was done

A full **276-MAG dereplicated combined reference** was built and the five
Nanopore samples were re-mapped to it, then the script was re-run against those
BAMs. Build script: `build_and_map_276ref.sh`.

- Reference: every `mag/drep/dereplicated_genomes/*.fa` header rewritten
  `>contig` → `>MAG_underscores::contig` (e.g. `CAN_5_bin_40::contig_7309_pilon`),
  concatenated → 14,633 contigs, 0 duplicate names.
- Aligner: minimap2 v2.30 `-ax map-ont --secondary=no` (primary alignments only).
- Sort/index/count: samtools v1.23.1.
- Reads: `processed/reads/nanopore/CAN_{1..5}.filt.fastq.gz` (filtlong-filtered ONT).
- Reference + BAMs were written to local NVMe `/scratch/afreiburger/nosz_allbins/`
  and are **not** committed (large). Re-run `build_and_map_276ref.sh` to regenerate.
- `compute_nosz_rpkm.py` `CONFIG.BAM_DIR` was set to `/scratch/afreiburger/nosz_allbins`.

Result: **naming check OK — all 136 loci across all 122 bins found**; every locus
has RPKM > 0 in ≥1 sample.

## Per-sample total primary-mapped reads (RPKM denominator)

| sample | full 276-ref (`--secondary=no`) | focal-ref note (README §5) |
|--------|---------------------------------:|---------------------------:|
| CAN_1 | 6,509,351  | ≈ 6.64M |
| CAN_2 | 5,224,260  | ≈ 5.05M |
| CAN_3 | 7,599,942  | ≈ 7.35M |
| CAN_4 | 7,156,876  | ≈ 9.26M |
| CAN_5 | 8,461,308  | ≈ 12.33M |

These differ from the focal-reference totals quoted in the README sanity check
because (a) the reference is the full 276-MAG catalogue rather than 21 focal MAGs,
and (b) this run counts **primary alignments only** (`-F 2308`), whereas the
focal `idxstats` totals included secondary/supplementary records. Each clade/bin
RPKM here uses its own sample's full-reference primary-mapped total as the
denominator, so the table is internally consistent; absolute RPKM values are
therefore **not** directly comparable to the original 13-MAG focal figure (the
all-community quantification is the more correct one for an all-bins comparison —
reads are assigned to their true source MAG rather than forced onto a focal
homolog).

## Outputs (`nosz_rpkm_out/`)

| file | content |
|------|---------|
| `nosz_locus_rpkm.tsv` | 136 nosZ loci × sample: reads, RPKM, depth-based copy number |
| `nosz_rpkm_per_sample_claded.tsv` | per bin × sample, long format, `ko_group ∈ {nosZ_cladeI, nosZ_cladeII}` (35 + 89 = 124 bin-rows; same layout as `gene_ab_figure/data/gene_rpkm_per_sample_claded.tsv`) |
| `nosz_copynumber_per_sample.tsv` | per bin × sample depth-based copy number, by clade |
| `sample_totals.tsv` | total primary-mapped reads per sample |

Clade tally over the 136 loci: 36 Clade I, 100 Clade II (no Clade III present in
the all-bins BED). Definitions unchanged from the README:
`RPKM = reads/(gene_kb)/(total_primary_mapped/1e6)`;
`copy_number = (reads/gene_len) / (bin_reads/bin_genome_len)`.
