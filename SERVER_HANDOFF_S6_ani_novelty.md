# Server handoff — Table S6 ANI novelty (+ Figure 1 prerequisites)

**Goal.** Finish the MAG novelty table (`Table S6_ANI` in `MAG_target_gene_filled.xlsx`)
and stage the inputs for the 62-MAG phylogenetic tree (Figure 1).

**What is already done on the workstation (no server needed):**
- `Table S6_ANI` is filled from `gtdbtk.bac120.summary.tsv` (already pulled): iterativeID,
  full GTDB taxonomy + 6 rank columns, closest reference (accession + species), ANI, AF,
  species-level call, novelty interpretation, and same-species co-member notes.
- Distribution: **91** "same species cluster", **1** edge ("same species, AF<0.60"),
  **184** "unresolved by ANI (GTDB species novel; `s__` unassigned)".

**Why the server is still needed.** GTDB-Tk's ANI screen only emits an ANI when the hit is
*already* within species range, so all 92 of its ANI values are ≥95 %. The 184 genomes with
no ANI hit are exactly the novel-species candidates — but GTDB-Tk gives them **no ANI number**,
so the rubric's *"putative novel species"* (ANI<95 & AF≥0.60) and *"divergent"* (ANI<95 & AF<0.60)
buckets cannot be populated, and the two coverage columns (O/P) are empty. **skani** vs the
GTDB r220 species reps fixes both. Parts A–C below.

---

## Conventions

```bash
# Pipeline roots (from the manuscript + scripts/build_mag_quality_table.py — verify with `ls`)
PROC=/scratch1/afreiburger/emilykin/processed
MAGDIR=$PROC/mag/drep/dereplicated_genomes          # the 276 final MAGs (*.fa)
GTDBTK=$PROC/mag/gtdbtk                              # classify/ align/ identify/
DREP=$PROC/mag/drep                                  # data_tables/ dereplicated_genomes/
OUT=$PROC/mag/skani_novelty;  mkdir -p $OUT

# GTDB-Tk reference data (r220). Find it:
echo "$GTDBTK_DATA_PATH"                             # if set, that's the release dir
#   else the manuscript says it is deployed at  refdata/gtdbtk/release220
GTDBREF=${GTDBTK_DATA_PATH:-/scratch1/afreiburger/emilykin/refdata/gtdbtk/release220}
```

Run long jobs in a **detached tmux** session (they survive logout):
```bash
tmux new -s s6ani        # … start the job …      Ctrl-b d  to detach;  tmux attach -t s6ani
```

Send results back to the workstation repo root (`/home/freiburger/Documents/EmilyKin/`), e.g.:
```bash
rsync -avP <files> freiburger@<workstation>:/home/freiburger/Documents/EmilyKin/
```

---

## Part A — Quick pulls (already computed; ~seconds)

These need no computation — just copy them back.

1. **dRep tables** — the receipt for the "same-species MAG-vs-MAG redundancy"
   (CAN_1_bin.174 vs CAN_5_bin.90). `Ndb.csv` holds the actual pairwise mutual ANI;
   `Cdb.csv` the cluster assignments; `Mdb.csv` the Mash primary buckets; `Wdb.csv` the winners.
   ```bash
   ls $DREP/data_tables/{Cdb,Ndb,Mdb,Wdb}.csv
   # send: $DREP/data_tables/Cdb.csv Ndb.csv Mdb.csv Wdb.csv
   ```

2. **120-marker MSA** — needed for the Figure 1 tree (Part C).
   ```bash
   ls $GTDBTK/align/gtdbtk.bac120.user_msa.fasta*       # .fasta or .fasta.gz
   # send: that user_msa file
   ```

3. **GTDB accession→taxonomy** — small; lets the workstation name skani's closest refs.
   ```bash
   ls $GTDBREF/taxonomy/gtdb_taxonomy.tsv               # 2 cols: RS_GCF_..\t d__..;..;s__..
   # send: gtdb_taxonomy.tsv   (rename to gtdb_taxonomy_r220.tsv on arrival)
   ```

---

## Part B — skani novelty run (recommended; completes the table)

Computes ANI + per-direction alignment fraction of all 276 MAGs to their nearest GTDB r220
species representative. This fills the "putative novel" / "divergent" calls for the 184
currently-unresolved MAGs and the Query/Reference coverage columns.

### B0. Get skani (pick one)
```bash
conda create -y -n skani -c bioconda -c conda-forge skani && conda activate skani
# — or a static binary —
# wget -qO skani https://github.com/bluenote-1577/skani/releases/latest/download/skani && chmod +x skani
skani -V
```

### B1. List the GTDB r220 reference genomes (robust to fastani/ vs skani/ layout)
```bash
# locate the reference genome FASTAs inside the GTDB-Tk data package
ls -d $GTDBREF/*/database 2>/dev/null
find $GTDBREF -maxdepth 4 -name '*_genomic.fna.gz' > $OUT/gtdb_r220_ref_paths.txt
wc -l $OUT/gtdb_r220_ref_paths.txt        # expect ~105k–113k reps (bac+arc, r220)
head -2 $OUT/gtdb_r220_ref_paths.txt
```
If that finds nothing, the genomes may be one level deeper — adjust `-maxdepth`, or use the
package's own manifest: `awk` the second column of `$GTDBREF/*/genome_paths.tsv` onto its dir.

### B2. Sketch the reference set once (reusable; ~20–40 min, 32 threads)
```bash
skani sketch -l $OUT/gtdb_r220_ref_paths.txt -o $OUT/gtdb_r220_sketch -t 32
```

### B3. Search the 276 MAGs against it (fast)
```bash
skani search -d $OUT/gtdb_r220_sketch -q $MAGDIR/*.fa \
    -o $OUT/mag_vs_gtdb_skani.tsv -t 32 --min-af 0 -n 3
wc -l $OUT/mag_vs_gtdb_skani.tsv          # >= 276 query rows (≤3 hits each)
head -3 $OUT/mag_vs_gtdb_skani.tsv        # cols: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
```
Notes: `--min-af 0` keeps the low-AF hits that define the "divergent" bucket; `-n 3` keeps the
top 3 per MAG (the workstation script takes the best-ANI hit). Align_fraction_* are **percent**.

### B4. Send back & fold into Table S6
```bash
# send: $OUT/mag_vs_gtdb_skani.tsv   (and gtdb_taxonomy_r220.tsv from Part A.3)
```
On the workstation:
```bash
cd /home/freiburger/Documents/EmilyKin
~/Documents/py_venv/bin/python scripts/build_table_s6_ani.py \
    --skani mag_vs_gtdb_skani.tsv --gtdb-tax gtdb_taxonomy_r220.tsv
# default: skani fills only the 184 no-ANI MAGs (GTDB-Tk stays authoritative for the 92).
# add --skani-all to recompute every row with skani for one-tool consistency.
```

---

## Part C — Figure 1 prerequisites (ON HOLD until the 62-MAG list arrives)

Figure: ML tree (RAxML, 100 bootstraps, 120 concatenated GTDB markers; Parks et al. 2018) of
the **62 selected MAGs + their closest reference genomes**, branches colored by taxonomy, red
circles = best ANI to closest reference (from Table S6), black circles on the 4 abundant/active
populations. **Blocked on:** the explicit list of 62 MAG bin names from the collaborator.

When the 62 list is in hand, the server needs:

1. **The user MSA** — `gtdbtk.bac120.user_msa.fasta` (Part A.2) gives the 120-marker
   concatenated alignment for the MAGs; subset it to the 62.
2. **Reference genomes on the tree.** Two options:
   - *Simplest, closest to the figure:* `gtdbtk de_novo_wf` seeded with the 62 MAGs and an
     outgroup, which co-aligns user genomes with GTDB reference genomes of the relevant taxa:
     ```bash
     gtdbtk de_novo_wf --genome_dir <dir_of_62_mags> --bacteria \
         --outgroup_taxon p__Patescibacteria \
         --out_dir $PROC/mag/tree62 --cpus 32 --extension fa
     ```
     (produces a decorated tree + the concatenated MSA of user + reference genomes).
   - *Curated/smaller:* keep only each MAG's closest reference(s) (from the Table S6
     `closest reference accession` + GTDB-Tk `other_related_references`), pull their bac120
     marker columns from `$GTDBREF/msa/*` (or the de_novo MSA), concatenate with the 62 user
     rows. The workstation can build this selection list once the 62 are fixed.
3. **RAxML, 100 rapid bootstraps, LG+GAMMA** on the concatenated protein MSA:
   ```bash
   raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -N 100 -m PROTGAMMALG \
       -s aln_62_plus_refs.faa -n fig1_62mags -T 16
   # → RAxML_bipartitions.fig1_62mags  (best tree with bootstrap support)
   ```
   (RAxML-NG equivalent: `raxml-ng --all --msa aln.faa --model LG+G --bs-trees 100`.)
4. **Send back:** `RAxML_bipartitions.fig1_62mags` (+ the MSA). The workstation renders the
   figure (taxonomy colors, ANI red circles from Table S6, black circles on the abundant set)
   with the existing figure tooling.

Proposed 4 black-circle "abundant & active" MAGs (confirm/adjust): *Ca.* Competibacter.139
(CAN_5_bin.106), *Ca.* Accumulibacter.8 (coasm_bin.185) & .23 (CAN_1_bin.98), Thauera.34
(CAN_4_bin.236).

---

## Send-back checklist

| File | From | Enables |
|------|------|---------|
| `Cdb.csv, Ndb.csv, Mdb.csv, Wdb.csv` | `$DREP/data_tables/` | confirms the same-species dedup explanation |
| `mag_vs_gtdb_skani.tsv` | Part B.3 | fills S6 "novel/divergent" calls + coverage cols O/P |
| `gtdb_taxonomy_r220.tsv` | `$GTDBREF/taxonomy/gtdb_taxonomy.tsv` | names skani's closest refs |
| `gtdbtk.bac120.user_msa.fasta` | `$GTDBTK/align/` | Figure 1 tree (when 62 list ready) |
| `RAxML_bipartitions.fig1_62mags` | Part C (later) | Figure 1 tree |
