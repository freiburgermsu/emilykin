#!/usr/bin/env python3
"""
Extract ppk1 ORFs (nucleotide) from our 7 'Ca. Accumulibacter' MAGs, following the
ppk1_Database workflow: pyrodigal gene calls -> pyhmmer hmmsearch vs ppk1.hmm
(top hit, E<=1e-50 first pass) -> pull the matching ORF's nucleotide CDS.

Writes ppk1_classify/our_ppk1.ffn and prints a per-MAG report.
"""
from pathlib import Path
import pyrodigal, pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from Bio.Seq import Seq

ROOT = Path("/home/freiburger/Documents/EmilyKin")
DREP = ROOT / "dereplicated_genomes"
HMM  = ROOT / "ppk1_classify/ppk1_Database/ppk1.hmm"
OUT  = ROOT / "ppk1_classify/our_ppk1.ffn"
MAGS = ["CAN_1_bin.98","CAN_4_bin.225","CAN_4_bin.64","coasm_bin.185",
        "coasm_bin.250","coasm_bin.347","coasm_bin.476"]

def read_fasta(p):
    seqs, cur, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if cur: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if cur: seqs[cur] = "".join(buf)
    return seqs

# ── call genes, collect proteins + nucleotide CDS ───────────────────────────
gf = pyrodigal.GeneFinder(meta=True)
amino = Alphabet.amino()
prot_index = {}          # protname -> (mag, nuc_seq, prot_seq, contig, length_aa)
digital = []
for mag in MAGS:
    contigs = read_fasta(DREP / f"{mag}.fa")
    ngenes = 0
    for contig, seq in contigs.items():
        genes = gf.find_genes(seq.encode())
        for i, g in enumerate(genes):
            ngenes += 1
            prot = g.translate().rstrip("*")
            nuc  = seq[g.begin-1:g.end]
            if g.strand == -1:
                nuc = str(Seq(nuc).reverse_complement())
            name = f"{mag}|{contig}|{i+1}"
            prot_index[name] = (mag, nuc, prot, contig, len(prot))
            digital.append(TextSequence(name=name.encode(), sequence=prot).digitize(amino))
    print(f"{mag}: {len(contigs)} contigs, {ngenes} ORFs")

# ── hmmsearch all proteins vs ppk1 HMM ──────────────────────────────────────
block = DigitalSequenceBlock(amino, digital)
with pyhmmer.plan7.HMMFile(HMM) as hf:
    hmm = hf.read()
hits_by_mag = {m: [] for m in MAGS}
for tophits in pyhmmer.hmmsearch([hmm], block, E=1e-10):
    for h in tophits:
        nm = h.name.decode() if isinstance(h.name, (bytes, bytearray)) else h.name
        mag = prot_index[nm][0]
        hits_by_mag[mag].append((h.score, h.evalue, nm))

# ── pick top hit per MAG, write nucleotide ORF ──────────────────────────────
print("\n=== ppk1 hits per MAG (sorted by bitscore) ===")
chosen = {}
for mag in MAGS:
    hl = sorted(hits_by_mag[mag], key=lambda t: -t[0])
    tag = "  (KO matrix said K00937=0)" if mag in ("CAN_1_bin.98","coasm_bin.476") else ""
    if not hl:
        print(f"{mag}: NO ppk1 hit (E<1e-10){tag}")
        continue
    print(f"{mag}: {len(hl)} hit(s){tag}")
    for sc, ev, nm in hl[:4]:
        L = prot_index[nm][4]
        print(f"    score={sc:7.1f}  E={ev:.1e}  len={L}aa  {nm}")
    sc, ev, nm = hl[0]
    pass1 = ev <= 1e-50
    chosen[mag] = (nm, sc, ev, pass1)

with open(OUT, "w") as f:
    for mag in MAGS:
        if mag in chosen:
            nm, sc, ev, _ = chosen[mag]
            _, nuc, prot, contig, L = prot_index[nm]
            f.write(f">{mag}_ppk1 len_nt={len(nuc)} len_aa={L} score={sc:.1f} E={ev:.1e} src={nm}\n")
            for j in range(0, len(nuc), 80):
                f.write(nuc[j:j+80] + "\n")
print(f"\nwrote {OUT}  ({len(chosen)} MAGs with ppk1)")
print("E<=1e-50 (README first-pass) pass:",
      [m for m in chosen if chosen[m][3]],
      "| weaker (1e-50..1e-10):", [m for m in chosen if not chosen[m][3]])
