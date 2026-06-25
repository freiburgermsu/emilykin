# Server-side processing — nosZ RPKM + depth-based copy number (all nosZ bins)

Everything except this step was done on the workstation (gene-calling, nosZ
identification by HMM, tree building, Clade I/II/III classification, and the
**genomic copy number** = count of nosZ loci per bin). The only thing that *must*
run on the server is read-counting from the per-sample BAM alignments — that's what
this folder is for.

This reproduces, for **all nosZ-bearing bins**, the same abundance metrics the
original figure had for the 13 selected MAGs:
- **RPKM** per nosZ locus and per bin, per sample, split by Clade I / II / III.
- **depth-based copy number** (gene depth ÷ genome depth) per bin, per sample, by clade.

---

## 0. Why the server

`compute_nosz_rpkm.py` needs the sample read alignments
`aligned_<sample>_sorted.bam` (reads mapped to the dereplicated-MAG combined
reference) and `samtools` — the exact inputs the workstation does **not** have but
that `gene_ab_figure/05_gene_copy_number.py` and `03_count_and_rpkm.py` already used
on poplar.

## 1. Files in this folder
| file | origin | role |
|------|--------|------|
| `compute_nosz_rpkm.py` | here | the script you run on the server |
| `allbins_nosz_loci_clade.bed` | workstation pipeline | every nosZ locus: BAM-style contig, coords, MAG, **Clade I/II/III** |
| `README_server_processing.md` | this file | directions |

`allbins_nosz_loci_clade.bed` columns (tab-sep):
`bam_contig  start0  end  locus_id  mag  clade  subtype  strand`
- `bam_contig` is already in the combined-reference / BAM naming
  `<mag_with_underscores>::<contig>` (e.g. `CAN_5_bin_40::contig_7309_pilon`).
- `start0` is 0-based (BED); the script converts to 1-based for samtools.

## 2. Prerequisites on the server
- `samtools` (any 1.x). The path used before was
  `/scratch1/afreiburger/emilykin/processed/.snakemake_envs/88fdb48d4d745c55ec2cd90b407de422_/bin/samtools`.
- The five sorted BAMs `aligned_CAN_1_sorted.bam … aligned_CAN_5_sorted.bam` in
  `/scratch1/afreiburger/emilykin/gene_ab_figure/data/` (the script will `samtools
  index` them if a `.bai` is missing).
- `python3` (standard library only — **no pip installs**).

## 3. Run it
```bash
# copy this folder to the server, then:
cd <this folder on the server>

# (a) open compute_nosz_rpkm.py and check the CONFIG block at the top:
#     SAMTOOLS, BAM_DIR, BAM_TMPL, SAMPLES  — defaults already match the paths above.

# (b) run
python3 compute_nosz_rpkm.py
```
Runtime: a few minutes (≈ n_loci × 5 `samtools view -c` calls).

## 4. What you get  (`nosz_rpkm_out/`)
| output | content |
|--------|---------|
| `nosz_locus_rpkm.tsv` | per nosZ locus × sample: read count, RPKM, copy number |
| `nosz_rpkm_per_sample_claded.tsv` | per bin × sample, long format, `ko_group ∈ {nosZ_cladeI, nosZ_cladeII, nosZ_cladeIII}` — **same layout as `gene_ab_figure/data/gene_rpkm_per_sample_claded.tsv`**, so it drops straight into that figure pipeline |
| `nosz_copynumber_per_sample.tsv` | per bin × sample depth-based copy number, by clade |
| `sample_totals.tsv` | total primary-mapped reads per sample (the RPKM denominator) |

### Definitions used
- `RPKM = reads_on_gene / (gene_len_kb) / (total_primary_mapped_reads / 1e6)`
- `copy_number = (reads_on_gene / gene_len) / (reads_on_bin_genome / bin_genome_len)`
- reads counted with `samtools view -c -F 2308` (primary mapped only: excludes
  unmapped + secondary + supplementary). Per-bin/genome totals from `samtools idxstats`.

## 5. Sanity checks the script prints
- per-sample total mapped reads (compare to your earlier run: CAN_1 ≈ 6.64M,
  CAN_2 ≈ 5.05M, CAN_3 ≈ 7.35M, CAN_4 ≈ 9.26M, CAN_5 ≈ 12.33M — from
  `05_classify_nosz.py`).
- a **naming check**: every BED contig must be found in the BAM. If it warns that
  contigs are missing, the BAM header naming differs from `combined_ref.fa`; tell me
  the BAM's header style (`samtools view -H <bam> | grep -m3 '^@SQ'`) and I'll
  re-emit the BED to match.

## 6. Bringing it back
Send `nosz_rpkm_out/` back and I'll merge it with the local clade + copy-number
table (`clade_classify/out/allbins_nosz_per_bin.tsv`) into one master nosZ table and,
if you want, regenerate the abundance/copy-number figures over all nosZ bins.

## 7. (Optional) same for ppk1
If you also want ppk1 RPKM for the Accumulibacter bins, say so — I'll emit a
`ppk1_loci_clade.bed` in the identical format and the **same script** computes it
(just point `BED=` at that file).
