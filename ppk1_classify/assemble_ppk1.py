#!/usr/bin/env python3
"""
(1) Relaxed re-check of the 2 MAGs with no ppk1 hit (E<1e-10): CAN_4_bin.225, coasm_bin.476.
(2) Assemble the combined nucleotide dataset for the tree:
       our ppk1 (+ any recovered partial) + DB refGenomes + clones + outgroup.
(3) Emit tip->clade metadata (refs from ppk1-database-info.csv; ours=QUERY; outgroup=OUTGROUP).
"""
from pathlib import Path
import csv
import pyrodigal, pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from Bio.Seq import Seq

ROOT = Path("/home/freiburger/Documents/EmilyKin")
PD   = ROOT / "ppk1_classify/ppk1_Database"
DREP = ROOT / "dereplicated_genomes"
OUTD = ROOT / "ppk1_classify"

def read_fasta(p):
    seqs, cur, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if cur: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if cur: seqs[cur] = "".join(buf)
    return seqs

# ── (1) relaxed re-check ────────────────────────────────────────────────────
gf = pyrodigal.GeneFinder(meta=True); amino = Alphabet.amino()
with pyhmmer.plan7.HMMFile(PD/"ppk1.hmm") as hf: hmm = hf.read()
recovered = {}
print("=== relaxed ppk1 search (E<=10) on the 2 missing MAGs ===")
for mag in ["CAN_4_bin.225","coasm_bin.476"]:
    contigs = read_fasta(DREP/f"{mag}.fa"); idx={}; dig=[]
    for contig,seq in contigs.items():
        for i,g in enumerate(gf.find_genes(seq.encode())):
            prot=g.translate().rstrip("*"); nuc=seq[g.begin-1:g.end]
            if g.strand==-1: nuc=str(Seq(nuc).reverse_complement())
            nm=f"{mag}|{contig}|{i+1}"; idx[nm]=(nuc,prot,len(prot)); dig.append(TextSequence(name=nm.encode(),sequence=prot).digitize(amino))
    hits=[]
    for th in pyhmmer.hmmsearch([hmm], DigitalSequenceBlock(amino,dig), E=10):
        for h in th:
            nm=h.name.decode() if isinstance(h.name,(bytes,bytearray)) else h.name
            hits.append((h.score,h.evalue,nm))
    hits.sort(key=lambda t:-t[0])
    if hits:
        for sc,ev,nm in hits[:3]:
            print(f"  {mag}: score={sc:6.1f} E={ev:.1e} len={idx[nm][2]}aa  {nm}")
        sc,ev,nm=hits[0]
        # include as partial only if a credible ppk1 fragment
        if sc>=80:
            recovered[mag]=(idx[nm][0], idx[nm][2], sc, ev)
            print(f"    -> RECOVER {mag} as partial ppk1 (score {sc:.0f})")
        else:
            print(f"    -> exclude {mag} (best score {sc:.0f} too low / not ppk1)")
    else:
        print(f"  {mag}: no hits even at E<=10  -> exclude")

# ── (2) assemble combined fasta ─────────────────────────────────────────────
ours = read_fasta(OUTD/"our_ppk1.ffn")
combined = OUTD/"combined_ppk1.ffn"
n_ref=n_clone=n_out=n_ours=0
with open(combined,"w") as out:
    for name,seq in ours.items():                    # our 5
        out.write(f">{name}\n{seq}\n"); n_ours+=1
    for mag,(nuc,L,sc,ev) in recovered.items():      # recovered partials
        out.write(f">{mag}_ppk1_partial\n{nuc}\n"); n_ours+=1
    for fn,tag in [("ppk1-refGenomes-coding-regions.fasta","ref"),
                   ("ppk1-clone-sequences.fasta","clone"),
                   ("outgroup-ppk1-coding-regions.fasta","out")]:
        for h,s in read_fasta(PD/"sequences"/fn).items():
            out.write(f">{h}\n{s}\n")
            if tag=="ref": n_ref+=1
            elif tag=="clone": n_clone+=1
            else: n_out+=1
print(f"\nwrote {combined}: ours={n_ours} refGenomes={n_ref} clones={n_clone} outgroup={n_out} "
      f"(total {n_ours+n_ref+n_clone+n_out})")

# ── (3) tip -> clade metadata ───────────────────────────────────────────────
loc2clade={}; loc2name={}
with open(PD/"ppk1-database-info.csv",encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        loc2clade[r["locus_tag"]]=r["Clade"]; loc2name[r["locus_tag"]]=r["name"]
meta=OUTD/"ppk1_tip_metadata.tsv"
with open(meta,"w") as m:
    m.write("tip\tcategory\tclade\tname\n")
    for name in ours: m.write(f"{name}\tQUERY\t?\t{name}\n")
    for mag in recovered: m.write(f"{mag}_ppk1_partial\tQUERY\t?\t{mag} (partial)\n")
    for fn,tag in [("ppk1-refGenomes-coding-regions.fasta","ref"),
                   ("ppk1-clone-sequences.fasta","clone"),
                   ("outgroup-ppk1-coding-regions.fasta","outgroup")]:
        for h in read_fasta(PD/"sequences"/fn):
            cl=loc2clade.get(h,"?"); nm=loc2name.get(h,h)
            cat = "OUTGROUP" if tag=="outgroup" else tag
            m.write(f"{h}\t{cat}\t{cl}\t{nm}\n")
print(f"wrote {meta}")
# unmatched check
unm=[h for fn in ["ppk1-refGenomes-coding-regions.fasta","ppk1-clone-sequences.fasta"]
     for h in read_fasta(PD/"sequences"/fn) if h not in loc2clade]
print(f"reference tips not found in CSV (clade=?): {len(unm)} {unm[:5]}")
