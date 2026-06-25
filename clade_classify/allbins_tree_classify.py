#!/usr/bin/env python3
"""
Build an all-bins nosZ protein tree (Chee+Orellana + He refs + every bin's nosZ
ORFs) and clan-classify each locus into Clade I/II/III — same marker-anchored,
rooting-robust method used for the 13-gene tree (partD), now for all nosZ bins.

In : out/allbins_nosz_query.faa , out/allbins_nosz_loci.tsv ,
     Chee_plus_Orellana_NosZ.prot.raw.faa , out/he_refs_degapped.faa
Out: out/allbins_nosz.aln , out/allbins_nosz.tree ,
     out/allbins_nosz_clades.tsv          (per locus: clade + evidence)
     out/allbins_nosz_loci_clade.bed      (for the server RPKM script)
     out/allbins_nosz_per_bin.tsv         (per bin: clade + copy number)
"""
from pathlib import Path
import csv, subprocess, os
from collections import Counter, defaultdict
from Bio import Phylo

CC = Path("/home/freiburger/Documents/EmilyKin/clade_classify"); OUT = CC / "out"
TOOLS = Path.home() / "Documents/emilykin_tools"   # stable location (survives session close)
MAFFT, FASTTREE = TOOLS / "mafft-linux64/mafft.bat", TOOLS / "FastTree"
I_ANCH = ["Pseudomonas","Paracoccus","Bradyrhizobium","Sinorhizobium","Rhodobacter",
          "Rhodopseudomonas","Shewanella","Marinobacter","Brucella","Azoarcus","Aromatoleum",
          "Thauera","Achromobacter","Ralstonia","Cupriavidus","Neisseria","Kingella","Roseobacter","Halomonas"]
II_ANCH = ["Wolinella","Anaeromyxobacter","Gemmatimonas","Chryseobacterium","Prevotella","Salinibacter",
           "Rhodothermus","Flavobacterium","Cytophaga","Dyadobacter","Haliscomenobacter","Gramella","Kordia",
           "Persicobacter","Chlorobaculum","Ignavibacterium","Bacteroides","Sphingobacterium"]

def read_fasta(p):
    s,cur,buf={},None,[]
    for l in open(p):
        if l.startswith(">"):
            if cur: s[cur]="".join(buf)
            cur=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if cur: s[cur]="".join(buf)
    return s

loci = {r["locus_id"]: r for r in csv.DictReader(open(OUT/"allbins_nosz_loci.tsv"), delimiter="\t")}
query_ids = set(loci)

# ── combined fasta (queries + refs) ─────────────────────────────────────────
comb = OUT/"allbins_nosz_combined.faa"
with open(comb,"w") as f:
    for k,v in read_fasta(OUT/"allbins_nosz_query.faa").items(): f.write(f">{k}\n{v}\n")
    # __r{i}/__h{i} suffix keeps names UNIQUE after truncation (FastTree rejects duplicates)
    for i,(k,v) in enumerate(read_fasta(CC/"Chee_plus_Orellana_NosZ.prot.raw.faa").items()):
        f.write(f">REF__{k.split()[0].replace('|','_')[:45]}__r{i}\n{v}\n")
    for i,(k,v) in enumerate(read_fasta(OUT/"he_refs_degapped.faa").items()):
        base=(k[4:] if k.startswith('HE__') else k)[:45]
        f.write(f">HE__{base}__h{i}\n{v}\n")

# ── align + tree ────────────────────────────────────────────────────────────
aln, tree_f = OUT/"allbins_nosz.aln", OUT/"allbins_nosz.tree"
env = dict(os.environ, MAFFT_BINARIES=str(TOOLS/"mafft-linux64/mafftdir/libexec"))
if aln.exists() and aln.stat().st_size > 0:
    print(f"reusing existing alignment {aln} (delete to force re-align)")
else:
    print("aligning (mafft --auto) ...")
    with open(aln,"w") as o: subprocess.run([str(MAFFT),"--auto","--quiet",str(comb)],stdout=o,check=True,env=env)
if tree_f.exists() and tree_f.stat().st_size > 0:
    print(f"reusing existing tree {tree_f} (delete to force re-build)")
else:
    print("building tree (FastTree -lg -gamma) ...")
    with open(aln) as i, open(tree_f,"w") as o, open(OUT/"allbins_nosz_fasttree.log","w") as e:
        subprocess.run([str(FASTTREE),"-lg","-gamma"],stdin=i,stdout=o,check=True,stderr=e)

# ── classify (clan-based) ───────────────────────────────────────────────────
tree = Phylo.read(tree_f,"newick")
is_ref=lambda n: n.startswith("REF__"); is_q=lambda n: n in query_ids
genus=lambda n: n.replace("REF__","").split("_")[0]
match=lambda n,ks: any(k.lower() in n.lower() for k in ks)
tips=[l.name for l in tree.get_terminals()]
I_tips=[t for t in tips if is_ref(t) and match(t,I_ANCH)]
II_tips=[t for t in tips if is_ref(t) and match(t,II_ANCH)]
# Root on a canonical Clade I anchor: midpoint mis-roots this larger tree (root lands
# inside Clade II, so its MRCA engulfs everything). Rooting on Pseudomonas/Paracoccus/
# Bradyrhizobium all give the SAME clean 330-tip Clade II clan (0 Clade I anchors) — a
# rooting-robust Clade I | Clade II bipartition.
og = (next((l for l in tree.get_terminals() if is_ref(l.name) and "pseudomonas_stutzeri" in l.name.lower()), None)
      or next((l for l in tree.get_terminals() if is_ref(l.name) and match(l.name, ["Pseudomonas","Paracoccus","Bradyrhizobium"])), None))
tree.root_with_outgroup(og)
mII=tree.common_ancestor([l for l in tree.get_terminals() if l.name in II_tips])
sII={l.name for l in mII.get_terminals()}
n_Iin=sum(t in set(I_tips) for t in sII)
print(f"Clade II clan: {len(sII)} tips | II-anchors {sum(t in set(II_tips) for t in sII)}/{len(II_tips)}"
      f" | I-anchors inside {n_Iin} | support {mII.confidence}")
assert n_Iin == 0, "Clade II clan not clean — check rooting anchor"
name2leaf={l.name:l for l in tree.get_terminals()}
def sister(q):
    for node in reversed([tree.root]+tree.get_path(name2leaf[q])[:-1]):
        refs=[t.name for t in node.get_terminals() if not is_q(t.name)]
        if refs:
            return node.confidence, [g for g,_ in Counter(genus(r) for r in refs if is_ref(r)).most_common(4)]
    return None,[]

clade_rows, bed_rows = [], []
for lid, r in loci.items():
    sub = r["subtype"]
    # tree placement is authority: inside Clade II clan -> II (even if the short GXHH motif
    # fired, e.g. truncated Bacteroidota fragments); only call III for loci OUTSIDE the
    # C-NosZ Clade II clan that also carry the L-NosZ (GXHH) CuZ motif.
    clade = "II" if lid in sII else ("III" if sub=="L-NosZ" else "I")
    supp, gen = sister(lid)
    clade_rows.append({"locus_id":lid,"mag":r["mag"],"clade":clade,"subtype":sub,
                       "support":"" if supp is None else f"{supp:.2f}","nearest_genera":";".join(gen)})
    bed_rows.append([r["bam_contig"], r["start0"], r["end"], lid, r["mag"], clade, sub, r["strand"]])

with open(OUT/"allbins_nosz_clades.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["locus_id","mag","clade","subtype","support","nearest_genera"],delimiter="\t")
    w.writeheader(); w.writerows(clade_rows)
with open(OUT/"allbins_nosz_loci_clade.bed","w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(bed_rows)

# ── per-bin merge with copy number ──────────────────────────────────────────
cn={r["MAG"]:r for r in csv.DictReader(open(OUT/"allbins_nosz_copynumber.tsv"),delimiter="\t")}
per=defaultdict(list)
for r in clade_rows: per[r["mag"]].append(r["clade"])
rows=[]
for mag in sorted(per):
    c=Counter(per[mag])
    rows.append({"MAG":mag,"nosZ_copies":len(per[mag]),"cladeI":c.get("I",0),"cladeII":c.get("II",0),
                 "cladeIII":c.get("III",0),"clade_call":";".join(sorted(set(per[mag]))),
                 "best_hmm_score":cn.get(mag,{}).get("best_score","")})
with open(OUT/"allbins_nosz_per_bin.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["MAG","nosZ_copies","cladeI","cladeII","cladeIII","clade_call","best_hmm_score"],delimiter="\t")
    w.writeheader(); w.writerows(rows)

nI=sum(r["cladeI"] for r in rows); nII=sum(r["cladeII"] for r in rows); nIII=sum(r["cladeIII"] for r in rows)
print(f"\n{len(rows)} nosZ bins | loci: Clade I={nI} Clade II={nII} Clade III={nIII}")
print("wrote allbins_nosz_clades.tsv, allbins_nosz_loci_clade.bed, allbins_nosz_per_bin.tsv")
# copy the server inputs next to the server script
import shutil
SP = CC/"server_processing"
shutil.copy(OUT/"allbins_nosz_loci_clade.bed", SP/"allbins_nosz_loci_clade.bed")
print(f"copied BED -> {SP}/allbins_nosz_loci_clade.bed")
