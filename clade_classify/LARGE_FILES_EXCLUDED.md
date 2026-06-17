# Large files excluded from git

Two Part B intermediates exceed GitHub's hard 100 MB per-file push limit, so they
are **kept on disk locally but excluded from version control** (see
`clade_classify/.gitignore`). Both are fully regenerable from inputs that *are*
tracked in this repo — nothing is lost.

| File | Size | Seqs | What it is |
|------|-----:|-----:|------------|
| `out/all_bins.fna` | 1,034,803,661 B (≈ 0.99 GiB) | 14,633 contigs | Nucleotide contigs of all 276 dereplicated MAGs, concatenated. The raw input for ORF calling. |
| `out/all_bins_orf.faa` | 343,758,000 B (≈ 0.33 GiB) | 935,110 ORFs | All predicted protein ORFs (pyrodigal *meta* mode) across the 276 MAGs. The searchable protein DB for the He 269NosZ HMM in Part B. |

## How to regenerate

Both derive from `dereplicated_genomes/` (the 276 MAG `*.fa` files), which **is**
tracked in this repo.

```bash
# 1. all_bins.fna — concatenate the dereplicated MAG contigs
cat dereplicated_genomes/*.fa > clade_classify/out/all_bins.fna

# 2. all_bins_orf.faa — call ORFs with pyrodigal (meta mode)
#    (edit the DREP/OUT paths at the top of the script for your machine)
python clade_classify/partB_step1_orf_call.py
```

`partB_step1_orf_call.py` reads each MAG, runs `pyrodigal.GeneFinder(meta=True)`,
and writes one FASTA record per ORF as `>{bin}::{contig}_{n}`. Output is
deterministic, so a regenerated `all_bins_orf.faa` matches the original.

## Note on the other excluded file

`clade_classify/He …等 - 2025 - A novel bacterial protein family … .pdf` (≈ 13 MB)
is the source paper (He et al. 2025, *Nature* 646:152). It is excluded via the
`*.pdf` rule (copyright; reading only, not a pipeline input).
