#!/usr/bin/env python3
"""
Classify our 5 'Ca. Accumulibacter' MAG ppk1 genes by placement in the
IQ-TREE ML tree (TVM+F+I+G4, SH-aLRT/UFBoot), rooted on the DB outgroup
(Dechloromonas aromatica + Rhodocyclus tenuis).

For each query: nearest ref-containing clan -> dominant DB ppk1 clade, the clan's
SH-aLRT/UFBoot support, and the nearest reference strain(s). Also checks whether the
DB clades themselves are monophyletic, for context.
"""
from pathlib import Path
import csv
from collections import Counter
from Bio import Phylo

OUTD = Path("/home/freiburger/Documents/EmilyKin/ppk1_classify")
tree = Phylo.read(OUTD/"ppk1_tree.treefile","newick")

# metadata
clade_of, name_of, cat_of = {}, {}, {}
for r in csv.DictReader(open(OUTD/"ppk1_tip_metadata.tsv"),delimiter="\t"):
    clade_of[r["tip"]]=r["clade"]; name_of[r["tip"]]=r["name"]; cat_of[r["tip"]]=r["category"]
# patch: Delftensis_ppk1 reference (clade IC per Petriglieri 2022)
for t in [l.name for l in tree.get_terminals()]:
    if t.lower().startswith("delftensis"): clade_of[t]="IC"; name_of[t]="Ca. A. delftensis"; cat_of[t]="ref"

QUERIES = [t for t in (l.name for l in tree.get_terminals()) if cat_of.get(t)=="QUERY"]
OUT_TIPS = [t for t in (l.name for l in tree.get_terminals()) if cat_of.get(t)=="OUTGROUP"]
REFS = {t:clade_of[t] for t in clade_of if cat_of.get(t) in ("ref","clone")}

# root on outgroup
tree.root_with_outgroup(tree.common_ancestor([l for l in tree.get_terminals() if l.name in OUT_TIPS]))

def sh_uf(node):
    lab = node.confidence if node.confidence is not None else node.name
    if lab is None: return ("","")
    s=str(lab)
    return tuple(s.split("/")) if "/" in s else (s,"")

# clade monophyly check
print("=== DB ppk1 clade monophyly in our tree ===")
for cl in sorted(set(REFS.values())):
    tips=[t for t,c in REFS.items() if c==cl]
    if len(tips)<2:
        print(f"  {cl:5s}: {len(tips)} tip (n/a)"); continue
    mrca=tree.common_ancestor([l for l in tree.get_terminals() if l.name in tips])
    clan=[l.name for l in mrca.get_terminals()]
    others=[clade_of.get(x,'?') for x in clan if x in REFS and REFS[x]!=cl]
    sh,uf=sh_uf(mrca)
    print(f"  {cl:5s}: {len(tips)} refs, clan={len(clan)} tips, intruders={len(others)} {Counter(others) if others else ''} SHaLRT/UFB={sh}/{uf}")

# place each query
name2leaf={l.name:l for l in tree.get_terminals()}
def place(q):
    path=tree.get_path(name2leaf[q]); anc=[tree.root]+path[:-1]
    for node in reversed(anc):
        refs=[t.name for t in node.get_terminals() if t.name in REFS]
        if refs:
            clades=Counter(REFS[r] for r in refs)
            dom,n=clades.most_common(1)[0]
            sh,uf=sh_uf(node)
            nearest=[(r,name_of.get(r,r),REFS[r]) for r in refs][:6]
            return dom, n/len(refs), clades, (sh,uf), nearest, len(refs)
    return "?",0,Counter(),("",""),[],0

print("\n=== QUERY ppk1 PLACEMENTS ===")
rows=[]
for q in sorted(QUERIES):
    dom,purity,clades,(sh,uf),nearest,nref=place(q)
    rows.append((q,dom,purity,clades,sh,uf,nearest,nref))
    print(f"\n{q}")
    print(f"  -> ppk1 clade = {dom}   (clan purity {purity:.0%}, {nref} refs in clan; SH-aLRT/UFBoot={sh}/{uf})")
    print(f"  clan clade mix: {dict(clades)}")
    print(f"  nearest refs:")
    for r,nm,cl in nearest: print(f"     [{cl:4s}] {nm[:55]}")

import json
json.dump([{"query":q,"clade":dom,"purity":purity,"clade_mix":dict(cl2),
            "sh_alrt":sh,"ufboot":uf,"nearest":[{ "tip":r,"name":nm,"clade":c} for r,nm,c in nr],"n_refs_in_clan":nref}
           for (q,dom,purity,cl2,sh,uf,nr,nref) in rows],
          open(OUTD/"partG_placements.json","w"), indent=1)
print(f"\nwrote {OUTD/'partG_placements.json'}")
