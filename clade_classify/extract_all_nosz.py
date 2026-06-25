#!/usr/bin/env python3
"""
Extract nosZ from ALL 276 dereplicated MAGs and record genomic COORDINATES so the
locus set can be (a) clade-typed by tree placement (local) and (b) read-counted for
RPKM on the server. pyrodigal ORFs -> pyhmmer hmmsearch vs He 269NosZ HMM (E<=1e-30)
-> CuA/CuZ motif subtype (C-NosZ clade I/II vs L-NosZ clade III).

Writes (to out/):
  allbins_nosz_query.faa / .ffn   protein / nucleotide per nosZ locus (locus_id headers)
  allbins_nosz_loci.tsv           locus_id, mag, bam_contig, raw_contig, start0, end, strand, subtype, hmm_score, evalue
  allbins_nosz_copynumber.tsv     per-bin genomic copy number (C-NosZ / L-NosZ loci)
bam_contig uses the combined_ref.fa / BAM naming  {mag_with_underscores}::{contig}.
"""
from pathlib import Path
import csv, re
from collections import defaultdict
import pyrodigal, pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from Bio.Seq import Seq

ROOT = Path("/home/freiburger/Documents/EmilyKin")
DREP = ROOT / "dereplicated_genomes"
HMM  = ROOT / "clade_classify/out/He_269NosZ.hmm"
OUT  = ROOT / "clade_classify/out"
EVAL = 1e-30
CUA  = re.compile(r"C.{2}FC.{3}H.EM")
CUZc = re.compile(r"D.HH"); CUZl = re.compile(r"G.HH")

def read_fasta(p):
    s, cur, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if cur: s[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if cur: s[cur] = "".join(buf)
    return s

def subtype(aa):
    # CuZ DXHH=C-NosZ (clade I/II), GXHH=L-NosZ (clade III). CuZ is N-terminal of CuA.
    m = CUA.search(aa); region = aa[:m.start()] if m else aa
    c, l = bool(CUZc.search(region)), bool(CUZl.search(region))
    if not c and not l:
        c, l = bool(CUZc.search(aa)), bool(CUZl.search(aa))
    return "L-NosZ" if (l and not c) else "C-NosZ"

amino = Alphabet.amino()
with pyhmmer.plan7.HMMFile(HMM) as hf: hmm = hf.read()
gf = pyrodigal.GeneFinder(meta=True)
bins = sorted(DREP.glob("*.fa"))
print(f"scanning {len(bins)} bins (He 269NosZ HMM, E<={EVAL:.0e}) ...")

loci, prot_faa, nuc_ffn = [], [], []
for i, fa in enumerate(bins):
    mag = fa.stem; safe = mag.replace(".", "_")
    contigs = read_fasta(fa)
    idx, dig = {}, []
    for contig, seq in contigs.items():
        for g in gf.find_genes(seq.encode()):
            aa = g.translate().rstrip("*"); nuc = seq[g.begin-1:g.end]
            if g.strand == -1: nuc = str(Seq(nuc).reverse_complement())
            nm = f"{contig}:{g.begin}:{g.end}"
            idx[nm] = (aa, nuc, contig, g.begin, g.end, g.strand)
            dig.append(TextSequence(name=nm.encode(), sequence=aa).digitize(amino))
    hits = []
    for th in pyhmmer.hmmsearch([hmm], DigitalSequenceBlock(amino, dig), E=EVAL):
        for h in th:
            nm = h.name.decode() if isinstance(h.name,(bytes,bytearray)) else h.name
            aa, nuc, contig, b, e, strand = idx[nm]
            hits.append((h.score, h.evalue, contig, b, e, strand, subtype(aa), aa, nuc))
    for sc, ev, contig, b, e, strand, sub, aa, nuc in sorted(hits, key=lambda t:-t[0]):
        lid = f"{safe}__{contig}__g{b}"
        loci.append({"locus_id":lid, "mag":mag, "bam_contig":f"{safe}::{contig}",
                     "raw_contig":contig, "start0":b-1, "end":e,
                     "strand":"+" if strand==1 else "-", "subtype":sub,
                     "hmm_score":f"{sc:.1f}", "evalue":f"{ev:.1e}"})
        prot_faa.append(f">{lid}\n{aa}\n"); nuc_ffn.append(f">{lid}\n{nuc}\n")
    if (i+1) % 60 == 0: print(f"  ...{i+1}/{len(bins)}")

(OUT/"allbins_nosz_query.faa").write_text("".join(prot_faa))
(OUT/"allbins_nosz_query.ffn").write_text("".join(nuc_ffn))
cols = ["locus_id","mag","bam_contig","raw_contig","start0","end","strand","subtype","hmm_score","evalue"]
with open(OUT/"allbins_nosz_loci.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(loci)

per = defaultdict(lambda: {"C":0,"L":0,"best":0.0})
for r in loci:
    per[r["mag"]]["C" if r["subtype"]=="C-NosZ" else "L"] += 1
    per[r["mag"]]["best"] = max(per[r["mag"]]["best"], float(r["hmm_score"]))
cn = [{"MAG":m,"nosZ_copies":d["C"]+d["L"],"cnosz_I_II":d["C"],"lnosz_cladeIII":d["L"],
       "best_score":f"{d['best']:.1f}"} for m,d in sorted(per.items())]
with open(OUT/"allbins_nosz_copynumber.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cn[0]), delimiter="\t"); w.writeheader(); w.writerows(cn)

ko = set()
with open(ROOT/"ko_copy_number_matrix.csv") as f:
    for r in csv.DictReader(f):
        if r.get("K00376","0") != "0": ko.add(r["MAG"])
cnb = set(per)
print(f"\nnosZ loci: {len(loci)}  in {len(per)} bins "
      f"(C-NosZ {sum(d['C'] for d in per.values())}, L-NosZ {sum(d['L'] for d in per.values())})")
print(f"KO K00376 bins={len(ko)}  HMM nosZ bins={len(cnb)}  overlap={len(ko&cnb)}  "
      f"KO-only={len(ko-cnb)}  HMM-only={len(cnb-ko)}")
print(f"wrote allbins_nosz_query.faa/.ffn, allbins_nosz_loci.tsv, allbins_nosz_copynumber.tsv")
