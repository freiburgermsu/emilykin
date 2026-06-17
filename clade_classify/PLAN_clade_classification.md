# Plan — rigorous nosZ clade I / II / III classification (run on the server)

Goal: replace the provisional clade I/II call for this system's nosZ genes with a
curated‑reference classification, and additionally screen for the **novel clade III
L‑NosZ** (lactonase‑type N₂O reductase; He et al. 2025, *Nature* 646:152) that standard
K00376 annotation misses. Heavy steps (GraftM, optional read search) are meant to run on
a compute server; the light steps and figure integration run anywhere.

This supersedes the DNA‑HMM clade I/II call in `gene_ab_figure/05_classify_nosz.py`
(kept for provenance) for the I‑vs‑II decision, and adds clade III, which that script
did not consider.

---

## 1. Background — the three‑group model and its diagnostics

| Group | Enzyme | KO | Reference here |
|------|--------|----|----------------|
| Clade I  | C‑NosZ (Cu, "typical" denitrifier) | **K00376** | Chee+Orellana faa |
| Clade II | C‑NosZ (Cu, "atypical") | **K00376** | Chee+Orellana faa |
| Clade III | **L‑NosZ** (lactonase‑type, novel; <35% aa id to C‑NosZ) | **not K00376** | He GraftM 269NosZ package |

Diagnostic motifs (He 2025, Box 1) — usable as an independent, deterministic check:
- **CuA** (all N₂ORs): `CX₂FCX₃HXEM`  → regex `C.{2}FC.{3}H.EM`
- **CuZ**: **`DXHH`** in C‑NosZ (clade I/II) vs **`GXHH`** in L‑NosZ (clade III) → `D.HH` vs `G.HH`
- Supporting CuZ motifs (both types): `XXHX`, `PHG`, `GPLH`, `XXGH`, `EPH`.

So: an N₂OR‑family hit with the CuA motif **and** a `GXHH` (not `DXHH`) CuZ site, at <35%
identity to canonical C‑NosZ, is L‑NosZ (clade III). Phylogenetic placement on the He tree
is the authority; the motif test is fast corroboration.

---

## 2. Inputs

In‑repo (`clade_classify/`):
- `Chee_plus_Orellana_NosZ.prot.raw.faa` — 336 C‑NosZ proteins (clade I + II; Chee clade I/II + Orellana typical/atypical). **Not** clade‑labeled per sequence — clade comes from the tree / marker taxa (see §5, §7).
- `graftM_NosZ.gpkg.refpkg.zip` — He 269NosZ **refpkg** (aligned fasta + rerooted tree `tree44qkgy_5.tre` + `..._taxonomy.csv` + `..._seqinfo.csv` + phylo model). NB: this is the *refpkg*, not a full runnable `.gpkg` — rebuild the gpkg in §4.
- (PDF of He 2025 is reference reading; gitignored, not needed to run.)

Elsewhere in repo / on server:
- `gene_ab_figure/data/target_genes.bed` — the 13 K00376 loci (one per MAG).
- `gene_ab_figure/data/nosz_clades.tsv` — prior provisional clade I/II calls (to compare).
- `gene_ab_figure/data/combined_ref.fa` — concatenated dereplicated‑MAG contigs (`{safe_mag}::{contig}` headers); source for the 13 nosZ CDS and for bin ORFs.
- `dereplicated_genomes/*.fa` — 276 bins (clade III bins‑mode search target).
- **(server only)** co‑assembly + per‑sample **contigs** — preferred unbinned target.
- **(server only)** `meta/longreads/*.fastq.gz` — 118 GB; raw‑read target (heavy, last resort).

---

## 3. Server environment

GraftM and its placement dependencies are bioconda‑tier (not pip/uv‑installable):

```bash
mamba create -n graftm -c bioconda -c conda-forge \
    graftm pplacer hmmer mafft diamond fasttree dendropy krona
conda activate graftm
# Python helpers (motif test, figure integration) — same env or a second one:
mamba install -c conda-forge -c bioconda pyrodigal pyhmmer biopython matplotlib numpy
```

Compute (Threadripper‑class node): **bins ≈ 5–15 min**, **contigs ≈ minutes–1 h**,
**raw reads ≈ ½–3 days** (translated search of ~150–200 Gbp dominates; long reads are a
poor fit for GraftM — prefer contigs).

---

## 4. Rebuild the runnable He gpkg from the refpkg

```bash
cd clade_classify
unzip -o graftM_NosZ.gpkg.refpkg.zip -d gpkg_x
# graftM create builds the search HMM from the alignment and assembles a .gpkg.
# (Verify flag names against your graftM version: `graftM create -h`.)
graftM create \
  --alignment    gpkg_x/graftM_269NosZ_deduplicated_aligned.fasta \
  --rerooted_tree gpkg_x/tree44qkgy_5.tre \
  --taxonomy     gpkg_x/graftM_269NosZ_taxonomy.csv \
  --output       He_269NosZ.gpkg
```
If `graftM create` rejects the prebuilt tree, let it build its own (`--sequences` =
de‑aligned fasta + `--taxonomy`), or build just the search HMM:
`hmmbuild He_269NosZ.hmm gpkg_x/graftM_269NosZ_deduplicated_aligned.fasta`.

---

## 5. Part A — reclassify the K00376 C‑NosZ genes as clade I vs II

The 13 K00376 genes are all C‑NosZ; refine I‑vs‑II against Chee+Orellana.

```bash
# 5.1 extract + translate the 13 nosZ CDS (CDS coords already in target_genes.bed)
#     reuse the extraction in 05_classify_nosz.py, then translate longest ORF:
python - <<'PY'   # writes nosz13.faa
import csv
from gene_ab_figure import ...   # or inline: read combined_ref.fa, slice [start:end], 6-frame translate longest ORF (Bio.Seq)
PY
```
Preferred (placement): make a Chee+Orellana gpkg and place the 13 proteins, **or** place
them on the He tree and read the clade region:
```bash
graftM create --sequences Chee_plus_Orellana_NosZ.prot.raw.faa --taxonomy <I_II_labels.csv> \
              --output CheeOrellana_CNosZ.gpkg          # needs a tip→cladeI/II label file
graftM graft --forward nosz13.faa --graftm_package CheeOrellana_CNosZ.gpkg --output partA/
```
Fallback (tree, no labels needed up front): align + tree, then read clade by marker taxa:
```bash
mafft --add nosz13.faa --reorder gpkg_x/graftM_269NosZ_deduplicated_aligned.fasta > A.aln
fasttree A.aln > A.tre   # then assign each query to the clade I vs II monophyletic group:
#   clade I markers : Pseudomonas, Paracoccus, Bradyrhizobium, Azoarcus, Thauera, Shewanella
#   clade II markers: Wolinella, Anaeromyxobacter, Bacteroidota, Gemmatimonadota, Accumulibacter_cladeIIA
```
Output → `out/partA_clade_I_II.tsv` (gene, mag, clade, method, support); diff against
`gene_ab_figure/data/nosz_clades.tsv` to see which provisional calls change.

---

## 6. Part B — detect novel clade III L‑NosZ ("does our system carry it?")

L‑NosZ is not K00376, so search **everything**, bins first.

```bash
# 6.1 BINS (fast, do first) — give GraftM nucleotide contigs; it ORF-calls internally.
#     Concatenate bins with bin-prefixed contig names so hits trace back to a MAG:
mkdir -p partB && : > partB/all_bins.fna
for f in ../dereplicated_genomes/*.fa; do
  b=$(basename "$f" .fa)
  sed "s/^>/>${b}::/" "$f" >> partB/all_bins.fna
done
graftM graft --forward partB/all_bins.fna --graftm_package He_269NosZ.gpkg \
             --output partB/bins --search_method hmmsearch --force
# -> partB/bins/all_bins/all_bins_read_tax.tsv  (placement/taxonomy per hit ORF)
#    partB/bins/all_bins/*_orf.fa               (hit ORF amino-acid seqs)
```

```bash
# 6.2 MOTIF DIAGNOSIS on the hit ORFs — the decisive C-NosZ vs L-NosZ test (deterministic).
python motif_diagnose.py partB/bins/all_bins/*_orf.fa > out/partB_motif.tsv
```
`motif_diagnose.py`:
```python
import re, sys
from Bio import SeqIO
CUA  = re.compile(r'C.{2}FC.{3}H.EM')   # CuA — all N2OR
CUZc = re.compile(r'D.HH')              # CuZ — C-NosZ (clade I/II)
CUZl = re.compile(r'G.HH')              # CuZ — L-NosZ (clade III)
print("orf\tCuA\tCuZ\tcall")
for r in SeqIO.parse(sys.argv[1], 'fasta'):
    s = str(r.seq); m = CUA.search(s)
    tail = s[m.end():] if m else s       # CuZ sits C-terminal of CuA
    c, l = bool(CUZc.search(tail)), bool(CUZl.search(tail))
    call = ('L-NosZ_cladeIII' if l and not c else
            'C-NosZ_cladeI_II' if c else 'ambiguous')
    print(f"{r.id}\t{'CuA' if m else 'noCuA'}\t{'GXHH' if l and not c else 'DXHH' if c else '?'}\t{call}")
```
Note: bare `D.HH`/`G.HH` are short — confirm each L‑NosZ candidate by (a) its GraftM
placement landing in the **clade III / L‑NosZ subtree** (§7), and (b) <35% identity to the
Chee+Orellana C‑NosZ set (`hmmsearch`/`phmmer` score gap, or a quick alignment). For
robustness, read the CuZ residue from a `hmmalign` of the ORF to the He alignment rather
than the raw regex.

```bash
# 6.3 (optional, server) UNBINNED fraction — prefer assembled contigs over raw reads.
graftM graft --forward <coassembly+per_sample_contigs>.fa \
             --graftm_package He_269NosZ.gpkg --output partB/contigs --search_method hmmsearch
# Raw-read last resort (½–3 days; diamond search is faster than hmmsearch on reads):
# for s in CAN_1 CAN_2 CAN_3 CAN_4 CAN_5; do
#   graftM graft --forward ../meta/longreads/${s}_nanopore.fastq.gz \
#                --graftm_package He_269NosZ.gpkg --output partB/reads_$s --search_method diamond
# done
```

---

## 7. Mapping GraftM output → clade I / II / III

GraftM's taxonomy is **organismal**, not clade‑numbered. To call a placement "clade III",
label the 259 reference tips by clade once and reuse:
- **Clade III (L‑NosZ)** = the deep‑branching L‑NosZ subtree (He Fig. 4a, light‑blue;
  Desulfobacterota / Nitrospinota / Chloroflexota / Bacillota L‑NosZ).
- **Clade I** = Pseudomonadota C‑NosZ (yellow); **Clade II** = the remaining C‑NosZ (pink).
Derive the tip→clade table from the rerooted tree topology (root at the C‑NosZ/L‑NosZ split;
the L‑NosZ side = III; within C‑NosZ, the Pseudomonadota‑typical clade = I, rest = II), or
from the He supplementary tip annotations. A query's clade = the labeled clade of the
reference subtree it places into.

---

## 8. Outputs

```
clade_classify/out/
  partA_clade_I_II.tsv     # 13 K00376 genes: clade I/II (curated), vs prior call
  partB_motif.tsv          # every N2OR hit ORF: CuA, CuZ (DXHH/GXHH), call
  partB_lnosz_hits.tsv     # confirmed L-NosZ: bin, ORF, placement clade, %id-to-CNosZ
  summary.md               # what changed for I/II; whether clade III is present, in which bins
```

---

## 9. Figure integration (back on the workstation)

1. **If Part A changes any clade I/II call** → update the clade labels feeding
   `gene_ab_figure/data/gene_rpkm_per_sample_claded.tsv` (regenerate via the Part‑A result,
   then re‑run `gene_ab_figure/04_plot_heatmap.py`). The figure already has separate
   `nosZ Clade I` / `nosZ Clade II` rows.
2. **If clade III L‑NosZ is found in a figure MAG** → add a `nosZ Clade III (L‑NosZ)` row
   under the `Denitrification (N₂O→N₂)` bracket. Because L‑NosZ is **not** K00376, its RPKM
   is not in the current matrices: extract the L‑NosZ ORF's CDS coordinates, add them to a
   BED, and recompute per‑sample RPKM exactly as in `05_classify_nosz.py` §6 (per‑sample
   totals already recovered: CAN_1 6.64M, CAN_2 5.05M, CAN_3 7.35M, CAN_4 9.26M, CAN_5 12.33M).
3. Update `gene_ab_figure/METHODS_nosz_clades.md` to cite Chee, Orellana 2014, and He 2025,
   and describe the GraftM placement + CuA/CuZ motif criteria.

---

## 10. Quick checklist

- [ ] `conda activate graftm`; rebuild `He_269NosZ.gpkg` (§4)
- [ ] Part A: `nosz13.faa` → placement/tree → `out/partA_clade_I_II.tsv` (§5)
- [ ] Part B bins: `graftM graft` on `all_bins.fna` → motif test → confirm L‑NosZ (§6)
- [ ] (opt) Part B contigs/reads on server (§6.3)
- [ ] tip→clade table for interpretation (§7)
- [ ] integrate into figure + methods (§9)
