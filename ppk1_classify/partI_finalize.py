#!/usr/bin/env python3
"""Finalize: add %identity corroboration to the classification table and write SUMMARY_ppk1.md."""
from pathlib import Path
import csv
from Bio import AlignIO

OUTD=Path("/home/freiburger/Documents/EmilyKin/ppk1_classify")
aln=AlignIO.read(OUTD/"ppk1_aligned.afa","fasta"); seq={r.id:str(r.seq) for r in aln}
clade_of,name_of,cat_of={},{},{}
for r in csv.DictReader(open(OUTD/"ppk1_tip_metadata.tsv"),delimiter="\t"):
    clade_of[r["tip"]]=r["clade"]; name_of[r["tip"]]=r["name"]; cat_of[r["tip"]]=r["category"]
for t in seq:
    if t.lower().startswith("delftensis"): clade_of[t]="IC"; name_of[t]="Ca. A. delftensis"; cat_of[t]="ref"
def pid(a,b):
    m=n=0
    for x,y in zip(a,b):
        if x=="-" or y=="-": continue
        n+=1; m+=x==y
    return 100*m/n if n else 0
refs=[t for t in seq if cat_of.get(t) in ("ref","clone")]
def best(q):
    s=sorted(((pid(seq[q],seq[r]),r) for r in refs),reverse=True)[0]
    return s[0], name_of.get(s[1],s[1]), clade_of.get(s[1],"?")

place={r["MAG"]:r for r in csv.DictReader(open(OUTD/"ppk1_classification.tsv"),delimiter="\t")}
NOTE={"CAN_1_bin.98":"Type I, divergent/novel (max id ~88%, no established clade)",
      "CAN_4_bin.64":"Type I, unclassified (≈UBA11064, 99% id; ~88% to named clades)",
      "coasm_bin.185":"clade IA","coasm_bin.250":"clade IC (delftensis vs meliphilus: closest delftensis)",
      "coasm_bin.347":"clade IIA (phosphatis lineage, UW1 type strain)"}
rows=[]
for mag,r in place.items():
    q=mag+"_ppk1"; p,nm,cl=best(q)
    r["top_pid"]=f"{p:.1f}% ({cl})"; r["note"]=NOTE.get(mag,"")
    rows.append(r)
# the 2 Accumulibacter MAGs with no recoverable ppk1 (listed for completeness -> 7 bins total)
NO_PPK1={
 "CAN_4_bin.225":"Accumulibacter, ~92% size but 99 contigs (N50 85 kb); no ppk1 ORF assembled - gene likely split across a contig break (KO K00937=1 = sub-HMM fragment, not a full gene)",
 "coasm_bin.476":"Accumulibacter, ~55% complete (2.2 Mb, 4 contigs); ppk1 not recovered in the assembly",
}
for mag,reason in NO_PPK1.items():
    rows.append({"MAG":mag,"ppk1_type":"no ppk1 recovered","proposed_species":"",
                 "SH-aLRT/UFBoot":"","top_pid":"","nearest_reference":"","note":reason})
cols=["MAG","ppk1_type","proposed_species","SH-aLRT/UFBoot","top_pid","nearest_reference","note"]
with open(OUTD/"ppk1_classification.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter="\t"); w.writeheader()
    for r in sorted(rows,key=lambda x:x["MAG"]): w.writerow({k:r.get(k,"") for k in cols})
print("final table written")
for r in sorted(rows,key=lambda x:x["MAG"]):
    print(f"  {r['MAG']:14s} {r['ppk1_type']:22s} {r['proposed_species']:30s} id={r['top_pid']}")

SUMMARY=f"""# ppk1-based classification of CAN 'Ca. Accumulibacter' MAGs

Mirrors the Petriglieri et al. 2022 (mBio) reevaluation / McDaniel et al. ppk1_Database
workflow: identify ppk1 by HMM, build a nucleotide ppk1 tree with the curated reference
set + outgroup, and read each MAG's ppk1 clade off the tree; corroborate by % identity.

## Inputs & method
- **Reference DB:** elizabethmcd/ppk1_Database (cloned) — `ppk1.hmm`, 40 genome + 41 clone
  reference ppk1 sequences with clade designations (`ppk1-database-info.csv`), and the
  bundled outgroup.
- **Outgroup:** `outgroup-ppk1-coding-regions.fasta` = *Dechloromonas aromatica* +
  *Rhodocyclus tenuis* (the McMahon-lab outgroup shipped with the DB). NB: Petriglieri 2022
  itself rooted on three *Azonexus* (formerly *Dechloromonas*) isolates — the same lineage,
  so the root is equivalent.
- **Our ppk1:** 7 GTDB-*Accumulibacter* MAGs screened. Genes called with **pyrodigal**
  (meta), ppk1 found by **pyhmmer hmmsearch** vs `ppk1.hmm` (top hit, E≤1e-50). 5 MAGs carry
  a detectable ppk1; the matching nucleotide ORF was extracted.
    - `CAN_1_bin.98` had ppk1 (E=3e-153) **despite KO matrix K00937=0** (KO false negative).
    - `CAN_4_bin.225` (KO K00937=1) and `coasm_bin.476` had **no** ppk1 hit even at E≤10
      (fragmented/incomplete MAGs) and were excluded.
    - The "2-copy" MAGs (185/250/347) resolve to exactly one true ppk1 (HMM rejects the paralog).
- **Alignment:** MAFFT v7.526 `--auto` (88 nucleotide seqs × 2315 cols).
- **Tree:** IQ-TREE 2.4.0, ModelFinder `-m MFP` → **TVM+F+I+G4**, 1000 UFBoot + 1000 SH-aLRT,
  rooted on the outgroup. (Paper used MAFFT --auto + IQ-TREE MFP → GTR+F+I+G4, 100 bootstraps —
  equivalent.)
- **Corroboration:** pairwise nucleotide % identity (the DB/ANI criterion: ~90–100% within
  clade, ~80% between clades).

## Results — ppk1 clade and proposed species (Petriglieri 2022 Table 1)

| MAG | ppk1 type | proposed species | SH-aLRT/UFBoot | top %id (clade) |
|-----|-----------|------------------|:--------------:|:---------------:|
| coasm_bin.185 | **IA**  | *Ca.* Accumulibacter **regalis**     | 99.5/100 | 100.0% (IA) |
| coasm_bin.250 | **IC**  | *Ca.* Accumulibacter **delftensis**  | 100/100  | 99.8% (IC) |
| coasm_bin.347 | **IIA** | *Ca.* Accumulibacter **phosphatis**  | 81.9/99  | 99.4% (IIA) |
| CAN_4_bin.64  | Type I (unclassified) | — (skipped, unclear) | 100/100 | 99.1% (≈UBA11064, un) |
| CAN_1_bin.98  | Type I (novel/divergent) | — (skipped, unclear) | — | ~88% (no clade) |
| CAN_4_bin.225 | no ppk1 detected | — | — | — |
| coasm_bin.476 | no ppk1 detected | — | — | — |

**Notes**
- IA→regalis, IIB→propinquus, IIG→affinis are 1:1; IC, IIA (and IIC/IID/IIF) hold multiple
  species, so the species was set from the nearest reference strain (delftensis vs meliphilus
  for IC; phosphatis/UW1 vs aalborgensis for IIA). coasm_bin.250 is 99.8% to *delftensis* and
  coasm_bin.347 is 99.4% to the UW1 *phosphatis* type strain.
- `CAN_1_bin.98` and `CAN_4_bin.64` are Type I but fall outside every named sub-clade (they
  group with the DB's own unclassified *Accumulibacter* sp. UBA11064); species skipped per
  the rule "skip if unclear." `CAN_1_bin.98` is the most divergent (≤88% to any reference) —
  a candidate novel Type I lineage.

## Files
- `ppk1_Database/` — cloned reference database
- `our_ppk1.ffn`, `combined_ppk1.ffn`, `ppk1_aligned.afa` — query/combined seqs + MAFFT alignment
- `ppk1_tree.treefile` / `.iqtree` / `.contree` — IQ-TREE ML tree + report
- `ppk1_tree_annotated.nwk` — tree with clade-labelled tips (open in FigTree/iTOL)
- `ppk1_reference_tree.png/.pdf` — rendered tree, queries highlighted
- `ppk1_classification.tsv` — this table; `partG_placements.json` — placement detail
- `extract_ppk1.py`, `assemble_ppk1.py`, `partG..partI*.py` — pipeline scripts
"""
(OUTD/"SUMMARY_ppk1.md").write_text(SUMMARY)
print(f"\nwrote {OUTD/'SUMMARY_ppk1.md'}")
