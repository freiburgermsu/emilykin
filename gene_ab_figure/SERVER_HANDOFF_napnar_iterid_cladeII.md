# Server handoff — add nar/nap rows, iterativeID x-axis, and Clade II sub-clades to the gene-abundance figure

**Run this on the server that originally built `gene_ab_figure/` (poplar,
`/scratch1/afreiburger/emilykin/...`).** A local session scoped and validated all
three changes, but the *canonical* pipeline lives here: the original 21-MAG
focal reference, the bakta GFF3 gene coordinates, the KOfam/eggNOG annotations,
and the original `minimap2`/`samtools` conda envs. Finishing on the server
reproduces the existing figure's RPKM **exactly** and adds the new rows on the
same footing.

---

## ⚠️ 0. The one thing you must get right — the reference scope

The committed figure's RPKM was produced by mapping the long reads to a reference
built from **only the 21 selected MAGs** (`02_build_and_align.sh` concatenates
`selected_mag_list.txt`), **not** the full 276-MAG community. This matters:

- The local session re-quantified against the full 276-MAG `combined_ref` and the
  existing values changed a lot — conserved denitrification genes (*nirS*, *nosZ*,
  *norB*, *nirB*) lost 70–87 % of their reads because, in a 276-genome reference,
  those reads map to homologous copies in the **non-selected** MAGs instead of the
  focal one. Example (CAN_2_bin.6 *nirS*, sample CAN_2): committed **824** reads,
  21-MAG ref **813** (reproduces), 276-MAG ref **52** (reads leaked away).
- Verified empirically: **21-MAG ref median count ratio vs committed = 0.989**;
  276-MAG ref = 0.13–0.99 (non-uniform). So **use the 21-MAG focal reference**
  (`combined_ref.fa` from `02_build_and_align.sh`). Do **not** substitute a
  full-community reference unless you intend to re-quantify everything (see §8).

Also reproduce the original **counting**: `samtools view -c -F 4 BAM contig:start+1-end`
(records overlapping the gene), with the per-sample denominator from
`samtools flagstat` "mapped". That is what `03_count_and_rpkm.py` already does.

---

## 1. The three tasks

| # | Task | Touches |
|---|------|---------|
| A | Add **nar** (narGHI) and **nap** (napAB) rows, max-RPKM over subunits (same policy as nirB/D and norC) | `01`, BED, `03`/aggregation, `04` |
| B | Label the x-axis (columns) by **MAG iterativeID** | `04_plot_heatmap.py` |
| C | Assign **Clade II sub-clades** (NosZREF 1577 A–H) to the nosZ genes | `07_cladeII_subclades.py` (new) |

The local session committed helper scripts to the repo
(`/scratch1/afreiburger/emilykin/EmilyKin/gene_ab_figure/`); `git pull` to get
them. **NB:** on the server the *working* dir is
`/scratch1/afreiburger/emilykin/gene_ab_figure` (the `WORK` path hard-coded in
`01`/`03`), which is **separate** from the repo checkout. Reconcile paths as you go.

---

## 2. Task A — identify the nar/nap genes

The broad 42-KO KOfam search already scored nar/nap (`kofamscan/hmmsearch.tblout`).
Add these 5 KOs to `TARGET_KOS` in `01_select_mags_and_genes.py` and re-run it, so
they flow into `data/target_genes.bed` with bakta coordinates:

```python
# 01_select_mags_and_genes.py
TARGET_KOS = {'K00362','K00363','K15864','K00368','K02305','K00376',
              'K04561',                       # norB (already added for the norC row)
              'K00370','K00371','K00374',     # nar: narG, narH, narI
              'K02567','K02568'}              # nap: napA, napB
```

Per-KO thresholds (already in `kofamscan/target_ko_list.txt`; `01` reads them):

| KO | gene | threshold | score_type |
|----|------|----------:|-----------|
| K00370 | narG | 820.57 | full |
| K00371 | narH | 332.73 | full |
| K00374 | narI | 189.37 | **domain** |
| K02567 | napA | 758.93 | full |
| K02568 | napB | 31.83 | **domain** |

### Cross-check — the 20 nar/nap genes the local session found (9 of 21 MAGs)

Coordinates below were recovered locally with pyrodigal and **validated 38/38
exact against existing bakta BED coordinates**, so your bakta GFF3 should give the
same. Verify your `target_genes.bed` nar/nap rows match these (contig is
`{safe_mag}::{orig_contig}`, BED 0-based half-open):

| MAG | gene | KO | locus_tag | contig | start | end |
|-----|------|----|-----------|--------|------:|----:|
| CAN_1_bin.98 | napB | K02568 | KLJHLB_02899 | contig_1207_pilon | 2436672 | 2437140 |
| CAN_1_bin.98 | napA | K02567 | KLJHLB_02902 | contig_1207_pilon | 2439141 | 2441667 |
| CAN_2_bin.125 | napB | K02568 | JODJNM_00633 | contig_18767_pilon | 1659 | 2214 |
| CAN_2_bin.125 | napA | K02567 | JODJNM_00634 | contig_18767_pilon | 2177 | 4718 |
| CAN_2_bin.6 | napA | K02567 | HMFKDB_00707 | contig_4248_pilon | 694194 | 696750 |
| CAN_2_bin.6 | napB | K02568 | HMFKDB_00710 | contig_4248_pilon | 698623 | 699085 |
| CAN_3_bin.203 | napA | K02567 | NGIEFB_01348 | contig_30743_pilon | 39317 | 41807 |
| CAN_3_bin.203 | napB | K02568 | NGIEFB_01351 | contig_30743_pilon | 43803 | 44241 |
| CAN_5_bin.106 | narG | K00370 | LLKMPD_00053 | contig_197_pilon | 32124 | 35625 |
| CAN_5_bin.106 | narH | K00371 | LLKMPD_00055 | contig_197_pilon | 35887 | 37441 |
| CAN_5_bin.40 | narI | K00374 | GNOLDL_00945 | contig_2512_pilon | 45987 | 46713 |
| CAN_5_bin.40 | narH | K00371 | GNOLDL_00947 | contig_2512_pilon | 47411 | 48935 |
| CAN_5_bin.40 | narG | K00370 | GNOLDL_00948 | contig_2512_pilon | 48959 | 52751 |
| CAN_5_bin.70 | napB | K02568 | MDDBPK_01030 | contig_3740_pilon | 544443 | 544890 |
| CAN_5_bin.70 | napA | K02567 | MDDBPK_01031 | contig_3740_pilon | 544956 | 547470 |
| coasm_bin.185 | napA | K02567 | FDDLBN_02790 | contig_3003_pilon | 164338 | 166864 |
| coasm_bin.185 | napB | K02568 | FDDLBN_02793 | contig_3003_pilon | 168862 | 169330 |
| coasm_bin.55 | narG | K00370 | FCMEHI_03857 | contig_60669_pilon | 330801 | 334566 |
| coasm_bin.55 | narH | K00371 | FCMEHI_03858 | contig_60669_pilon | 334567 | 336121 |
| coasm_bin.55 | narI | K00374 | FCMEHI_03861 | contig_60669_pilon | 336963 | 337674 |

(nap = napAB in 6 MAGs: CAN_1_bin.98, CAN_2_bin.125, CAN_2_bin.6, CAN_3_bin.203,
CAN_5_bin.70, coasm_bin.185. nar = narGHI in 3: CAN_5_bin.106 [narGH only],
CAN_5_bin.40, coasm_bin.55. Genes are in operon order.)

**The committed `data/target_genes.bed` already contains these 20 nar/nap rows**
(merged by the local session, pyrodigal coords validated 38/38 against existing
bakta BED coords) — use it directly. Do **not** append `data/napnar_genes.bed`
again; it is the same 20 rows kept as a standalone record. If you prefer canonical
bakta coordinates, re-run `01` with the 5 KOs added (it regenerates
`target_genes.bed` from the GFF3, overwriting the merged file).

---

## 3. Task A — counting + aggregation (the figure's data file)

The figure reads `data/gene_rpkm_per_sample_claded.tsv`. Regenerate it so the
**existing rows are byte-for-byte reproduced** and two new rows are added.

**Counting:** if the per-sample BAMs from the original run still exist
(`data/aligned_CAN_{1..5}_sorted.bam` + `.bai`, 21-MAG ref), just count the new
nar/nap intervals over them — **no re-alignment needed**. Otherwise re-run
`02_build_and_align.sh` first (rebuilds the 21-MAG `combined_ref.fa` and the BAMs).

**Aggregation — the 8 figure rows and their KO sets (max RPKM over the row's genes
per MAG/sample):**

```
nar           max(K00370 narG, K00371 narH, K00374 narI)     # NEW
nap           max(K02567 napA, K02568 napB)                  # NEW
nirB/D        max(K00362, K00363)
nirS          K15864
nirK          K00368
norC          max(K02305 norC, K04561 norB)                  # see gotcha below
nosZ_cladeI   K00376 where nosz_clades.tsv clade == I
nosZ_cladeII  K00376 where nosz_clades.tsv clade == II
```

> **⚠️ Stale-script gotcha.** The committed `05_classify_nosz.py` §6 builds the
> norC row from **K02305 only** (`in_row(... 'norC' ... ko=='K02305')`), but the
> *current* figure's norC row is **max(K02305, K04561)** — norB was added later
> (10 MAGs, not 4). Do **not** regress this. Likewise `03_aggregate_rpkm_local.py`
> / `03_count_and_rpkm.py` `KO_GROUPS` are the old 4-group set. Use the row/KO map
> above. The repo's `gene_ab_figure/03c_aggregate_claded.py` already encodes
> exactly this 8-row policy and the nosZ clade split — adapt its count source to
> your BAM counts (it reads `counts_<sample>.json` = `{"total_mapped": N,
> "gene_counts": {gene_id: count}}`; easiest path: emit that JSON from
> `samtools view -c` + `flagstat` and run `03c` unchanged).

After regenerating, **verify** the existing rows match the committed file (they
must, on the 21-MAG ref):

```bash
git show HEAD:gene_ab_figure/data/gene_rpkm_per_sample_claded.tsv > /tmp/old.tsv
# diff the nirB/D, nirS, nirK, norC, nosZ_cladeI, nosZ_cladeII rows vs the new file;
# max abs diff should be ~0 (rounding). If not, the reference scope or counting differs.
```

---

## 4. Task A — plot rows (`04_plot_heatmap.py`)

Add a new **top** row-group. Prepend to `ROW_GROUPS`:

```python
ROW_GROUPS = [
    ('Nitrate reduction\n(NO₃→NO₂)',
     [('nar', 'nar  –  membrane-bound nitrate reductase (narGHI)'),
      ('nap', 'nap  –  periplasmic nitrate reductase (napAB)')]),
    ('DNRA\n(NO₂→NH₄)',
     [('nirB/D', 'nirB/D  –  NADH nitrite reductase')]),
    # ... existing groups unchanged ...
]
```

(The repo's `04_plot_heatmap.py` already has this edit plus Task B — `git pull` and
diff rather than re-editing by hand.)

---

## 5. Task B — iterativeID x-axis (`04_plot_heatmap.py`)

Label each MAG column by `mag_iterativeID` (from `data/taxonomy_labels.tsv`)
instead of `genus (short_id)`. Add a helper and use it for the bottom label:

```python
def iterid_of(m):
    it = meta.get(m, {}).get('iterid', '')
    return it.replace('Ca_', 'Ca. ') if it else short_id(m)
```
```python
# in the per-MAG bottom-label loop, replace the genus/short_id text with:
    ax.text(bx + n_s * SUB_W / 2, total_y + 0.5,
            iterid_of(m), ha='right', va='top',
            rotation=45, rotation_mode='anchor', fontsize=6.8,
            fontstyle='italic', clip_on=False)
```

`meta[m]['iterid']` is already loaded from `taxonomy_labels.tsv`
(`mag_iterativeID`). Examples: `Ca. Accumulibacter.23`, `Dokdonella.14`,
`UBA5066_m.1`, `Denitratisoma.2`. Update the caption to say columns are labelled by
MAG iterativeID and that nar/nap/nirB/D/norC cells are the max over their subunit
genes.

---

## 6. Task C — Clade II sub-clades

Run `gene_ab_figure/07_cladeII_subclades.py` (path-portable; `git pull` to get it).
It nhmmer-scores each Clade-II nosZ gene against the 8 NosZREF DNA HMMs (same
machinery as `05_classify_nosz.py`) and writes the best sub-clade + runner-up +
margin to `data/nosz_cladeII_subclades.tsv`. It needs `data/combined_ref.fa`,
`data/target_genes.bed`, `data/nosz_clades.tsv`, and `literature_claudeI_II/`.

**Expected result (local; should reproduce):**

| sub-clade | n | MAGs (genus) | confidence |
|-----------|--:|--------------|-----------|
| II.D | 7 | Accumulibacter ×2, Desulfobacillus/Denitratisoma, Azonexus, JAEULV01, JAAEKA01, Leptovillus | high (margins 85–160); coasm_bin.260 & .481 medium (~12) |
| II.G | 2 | Giesbergeria (high, 60), UBA5066 (medium, 14) | — |
| II.C | 2 | JJ008, SpSt-398 (Terrimonas) | **low** — II.C vs II.B margin 0.8 |
| II.B | 1 | OLB5 | **low** — margin 2.2 |

Caveat: II.B and II.C are both Bacteroidetes; the three Bacteroidetes genes sit at
the II.B/II.C boundary (<2.5 bits) → report as "Clade II Bacteroidetes (II.B/C)".
II.D (Dechloromonas-type Betaproteobacterial denitrifiers) and II.G are robust.
See `METHODS_nosz_clades.md` §6 (added by the local session).

---

## 7. Recommended execution order on the server

```
1.  git pull                              # get 04 (rows+iterid), 03c, 07, napnar_genes.bed, METHODS §6
2.  add the 5 KOs to 01 → re-run 01       # OR append data/napnar_genes.bed to target_genes.bed
3.  confirm reference = 21 selected MAGs  # combined_ref.fa from 02_build_and_align.sh
4.  counts:
      if data/aligned_CAN_*_sorted.bam exist:  count target_genes.bed intervals with
            samtools view -c -F 4 ; total = samtools flagstat 'mapped'
      else:  run 02_build_and_align.sh  (re-aligns; 21-MAG ref)
5.  aggregate → gene_rpkm_per_sample_claded.tsv   # 8 rows, max policy, norC=max(K02305,K04561)
6.  VERIFY existing 6 rows == committed (max abs diff ~0)     # §3
7.  run 04_plot_heatmap.py → gene_abundance_figure.png/.pdf
8.  run 07_cladeII_subclades.py → nosz_cladeII_subclades.tsv
9.  eyeball the figure (nar/nap top group present; iterativeID labels; existing rows unchanged)
10. commit + push   (guard: no file >100 MB — combined_ref.*, BAMs, *.mmi stay gitignored)
```

---

## 8. If you (intentionally) want the full-community quantification instead

Mapping to all 276 dereplicated MAGs is arguably *more accurate* (a read is
assigned to its true source genome rather than force-mapped onto one of the 21
focal MAGs), and it is what `METHODS_nosz_clades.md` / the `03_aggregate_*`
docstrings actually describe. But it **changes the existing figure** (conserved
genes drop sharply, as in §0). If that re-quantification is desired, do it
deliberately and as a separate, documented change to all rows — not as a
side effect of adding nar/nap. Default for this task: **stay on the 21-MAG ref**
so the only changes are the two new rows, the x-axis labels, and the sub-clade table.

---

## 9. Files the local session prepared (in the repo, `git pull`)

| file | purpose |
|------|---------|
| `04_plot_heatmap.py` | already edited: nar/nap top group + iterativeID x-axis + caption |
| `03c_aggregate_claded.py` | clean 8-row aggregation (max policy, norC=max(K02305,K04561), nosZ clade split) from per-sample counts JSON |
| `07_cladeII_subclades.py` | Clade II sub-clade assignment (nhmmer vs NosZREF 1577 A–H) |
| `01b_napnar_coords.py` | local coord recovery + cross-check (server should prefer bakta GFF3) |
| `count_stream_sam.py` | SAM-streaming read counter (only if you stream minimap2 instead of samtools-counting BAMs) |
| `data/napnar_genes.bed` | the 20 nar/nap BED rows (coords above) |
| `data/nosz_cladeII_subclades.tsv` | Clade II sub-clade result |
| `METHODS_nosz_clades.md` | §6 added: sub-clade method + result |
| `08_gene_copy_per_mag.py` | genomic gene-copy inventory — copies of each KO per MAG, nosZ split Clade I/II; reference-independent (counts BED loci), no server action needed |
| `data/gene_copy_per_mag.tsv` / `.md` | the copy-inventory table (already complete) |

**Do not** trust the local `data/counts_*.json`, `data/gene_rpkm_per_sample_claded.tsv`,
or `data/read_counts_raw.tsv` if present in the working tree — the local session
generated those against the **276-MAG** ref for diagnosis; regenerate them on the
21-MAG ref here.
