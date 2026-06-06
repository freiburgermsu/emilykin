#!/usr/bin/env python
"""
prep_inputs.py — turn the QIIME2 artifacts into the inputs the embedding-hit pipeline expects.

Reads (from ../qiime, the merged-DADA2 set):
  rep_seqs_merged_dada2.qza  -> data/dna-sequences.fasta   (the ASV sequences)
  taxonomy_merged.qza        -> data/taxonomy.tsv          (MiDAS k__..;g__..;s__.. strings)
  table_merged_dada2.qza     -> data/feature-table.biom    (per-ASV per-sample counts, HDF5)

Writes (to ./inputs):
  asvs.fasta      copy of the rep-seqs
  taxonomy.csv    seq,Kingdom,Phylum,Class,Order,Family,Genus,Species,rel_ab   (rel_ab = max % across samples)
  asv_IDs.csv     ASV, MiDAS Taxonomy                                          (joined lineage string)

A .qza is just a zip: <uuid>/data/<file>. No qiime2 install required.
"""
from __future__ import annotations
import zipfile, csv
from pathlib import Path
import numpy as np, scipy.sparse as sp, h5py

QIIME = Path(__file__).resolve().parent.parent / "qiime"
OUT = Path(__file__).resolve().parent / "inputs"
RANK_LETTER = {"k": "Kingdom", "d": "Kingdom", "p": "Phylum", "c": "Class",
               "o": "Order", "f": "Family", "g": "Genus", "s": "Species"}
RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def read_member(qza, suffix):
    z = zipfile.ZipFile(QIIME / qza)
    name = [n for n in z.namelist() if n.endswith(suffix)][0]
    return z.read(name)


def parse_taxon(s):
    """'k__Bacteria; p__Proteobacteria; ...; g__Thauera; s__' -> {rank: name}."""
    out = {r: "" for r in RANKS}
    for tok in s.split(";"):
        tok = tok.strip()
        if len(tok) >= 3 and tok[1:3] == "__":
            rank = RANK_LETTER.get(tok[0].lower())
            if rank:
                out[rank] = tok[3:].strip()
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # 1. rep-seqs fasta
    (OUT / "asvs.fasta").write_bytes(read_member("rep_seqs_merged_dada2.qza", "dna-sequences.fasta"))

    # 2. per-ASV abundance from the BIOM feature table
    biom = OUT / "feature-table.biom"
    biom.write_bytes(read_member("table_merged_dada2.qza", "feature-table.biom"))
    with h5py.File(biom, "r") as h:
        obs = [x.decode() for x in h["observation/ids"][:]]
        samp = [x.decode() for x in h["sample/ids"][:]]
        M = sp.csc_matrix((h["sample/matrix/data"][:], h["sample/matrix/indices"][:],
                           h["sample/matrix/indptr"][:]), shape=(len(obs), len(samp)))   # obs x samp
    samp_tot = np.asarray(M.sum(0)).ravel()
    rel = M.multiply(1.0 / np.maximum(samp_tot, 1))
    max_relab = {obs[i]: float(np.asarray(rel.getrow(i).todense()).ravel().max() * 100) for i in range(len(obs))}

    # 3. taxonomy.tsv -> ranks
    tsv = read_member("taxonomy_merged.qza", "taxonomy.tsv").decode().splitlines()
    rows, idrows = [], []
    for line in tsv[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        asv, taxon = parts[0], parts[1]
        r = parse_taxon(taxon)
        rows.append({"seq": asv, **r, "rel_ab": round(max_relab.get(asv, 0.0), 6)})
        idrows.append({"ASV": asv, "MiDAS Taxonomy": " ".join(r[k] for k in RANKS if r[k])})

    with open(OUT / "taxonomy.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["seq"] + RANKS + ["rel_ab"])
        w.writeheader(); w.writerows(rows)
    with open(OUT / "asv_IDs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ASV", "MiDAS Taxonomy"])
        w.writeheader(); w.writerows(idrows)

    g_assigned = sum(1 for r in rows if r["Genus"] and not r["Genus"].lower().startswith("midas"))
    print(f"ASVs={len(rows)} | samples={len(samp)} | genus(real-named)={g_assigned} "
          f"({100*g_assigned/len(rows):.0f}%) | wrote inputs/asvs.fasta, taxonomy.csv, asv_IDs.csv")


if __name__ == "__main__":
    main()
