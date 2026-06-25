#!/usr/bin/env python3
"""
Emit ppk1 deliverables:
  ppk1_tree_annotated.nwk : tree with tips relabelled by clade; queries marked.
  ppk1_reference_tree.png/.pdf : full 88-tip ML tree, queries in red, refs by Type.
  ppk1_classification.tsv : our 5 MAGs -> ppk1 clade -> proposed species (Petriglieri 2022).
"""
from pathlib import Path
import csv, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

OUTD = Path("/home/freiburger/Documents/EmilyKin/ppk1_classify")
TYPE_I  = {"IA","IB","IC","ID","IE"}
TYPE_II = {"IIA","IIB","IIC","IID","IIE","IIF","IIG","IIH","II-I"}
# clade -> proposed species (Petriglieri 2022 Table 1); [] = none proposed
CLADE_SPECIES = {"IA":["regalis"],"IB":["appositus","adiacens"],"IC":["meliphilus","delftensis"],
 "IIA":["aalborgensis","phosphatis"],"IIB":["propinquus"],"IIC":["contiguus","vicinus","cognatus"],
 "IID":["proximus","necessarius"],"IIF":["iunctus","adjunctus","similis","conexus"],"IIG":["affinis"]}
# strain token -> species (to resolve multi-species clades by nearest reference)
STRAIN_SPECIES = {"UW3":"regalis","BA-93":"regalis","CANDO_1":"regalis","UW4":"regalis","UW8":"regalis",
 "BA-92":"appositus","HKU1":"adiacens","UW-LDO":"meliphilus","delftensis":"delftensis","SBR_S":"delftensis",
 "UW1":"phosphatis","UW9":"phosphatis","UW5":"phosphatis","aalborgensis":"aalborgensis",
 "SBR_L":"contiguus","SK-01":"vicinus","UBA5574":"vicinus","Bin19":"cognatus","HKU2":"cognatus",
 "SK-02":"cognatus","SSA1":"cognatus","UW11":"cognatus","UW6":"cognatus","UBA2327":"iunctus",
 "SK-12":"adjunctus","SCELSE-1":"similis","SSB1":"similis","UW13":"conexus","UW7":"conexus",
 "Fred_BAT3C.720":"affinis","EsbW_BATAC.285":"proximus","UW12":"necessarius"}

tree = Phylo.read(OUTD/"ppk1_tree.treefile","newick")
clade_of,name_of,cat_of={},{},{}
for r in csv.DictReader(open(OUTD/"ppk1_tip_metadata.tsv"),delimiter="\t"):
    clade_of[r["tip"]]=r["clade"]; name_of[r["tip"]]=r["name"]; cat_of[r["tip"]]=r["category"]
for t in [l.name for l in tree.get_terminals()]:
    if t.lower().startswith("delftensis"): clade_of[t]="IC"; name_of[t]="Ca. A. delftensis"; cat_of[t]="ref"
out=[l for l in tree.get_terminals() if cat_of.get(l.name)=="OUTGROUP"]
tree.root_with_outgroup(tree.common_ancestor(out))

place={p["query"]:p for p in json.load(open(OUTD/"partG_placements.json"))}
# refine: the two 'un' queries are Type I (unclassified) per neighbourhood analysis
REFINE={"CAN_1_bin.98_ppk1":"Type I (unclassified)","CAN_4_bin.64_ppk1":"Type I (unclassified)"}

def species_for(clade, nearest):
    if clade not in CLADE_SPECIES: return ""
    sp=CLADE_SPECIES[clade]
    if len(sp)==1: return sp[0]
    for n in nearest:                       # multi-species clade: use nearest strain
        for tok,s in STRAIN_SPECIES.items():
            if tok.lower() in n["name"].lower(): return s
    return ""                                # unresolved within clade

# ── classification table ────────────────────────────────────────────────────
rows=[]
for q,p in sorted(place.items()):
    clade=p["clade"]; disp=REFINE.get(q,clade)
    sp = "" if q in REFINE else species_for(clade,p["nearest"])
    rows.append({"MAG":q.replace("_ppk1",""),"ppk1_type":disp,
                 "proposed_species": (f"Ca. Accumulibacter {sp}" if sp else ""),
                 "SH-aLRT/UFBoot":f"{p['sh_alrt']}/{p['ufboot']}",
                 "nearest_reference":"; ".join(n["name"] for n in p["nearest"][:3])})
with open(OUTD/"ppk1_classification.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["MAG","ppk1_type","proposed_species","SH-aLRT/UFBoot","nearest_reference"],delimiter="\t")
    w.writeheader(); w.writerows(rows)
print("=== ppk1_classification.tsv ===")
for r in rows: print(f"  {r['MAG']:14s} {r['ppk1_type']:22s} {r['proposed_species']:30s} {r['SH-aLRT/UFBoot']}")

# ── annotated newick ────────────────────────────────────────────────────────
t2=Phylo.read(OUTD/"ppk1_tree.treefile","newick"); t2.root_with_outgroup(t2.common_ancestor([l for l in t2.get_terminals() if l.name in [o.name for o in out]]))
for l in t2.get_terminals():
    c=clade_of.get(l.name,"?");
    if cat_of.get(l.name)=="QUERY":
        l.name=f"[Q:{REFINE.get(l.name,place.get(l.name,{}).get('clade','?'))}]{l.name}"
    elif cat_of.get(l.name)=="OUTGROUP": l.name=f"[OUT]{name_of.get(l.name,l.name)[:30]}"
    else: l.name=f"[{c}]{name_of.get(l.name,l.name)[:35]}"
Phylo.write(t2,OUTD/"ppk1_tree_annotated.nwk","newick")
print(f"\nwrote {OUTD/'ppk1_tree_annotated.nwk'}")

# ── figure ──────────────────────────────────────────────────────────────────
def tcol(n):
    if cat_of.get(n)=="QUERY": return "#12355B"
    if cat_of.get(n)=="OUTGROUP": return "black"
    c=clade_of.get(n,"?")
    return "#1f77b4" if c in TYPE_I else "#2ca02c" if c in TYPE_II else "#ff7f0e"
def lab(cl):
    n=cl.name
    if not n: return None
    c=clade_of.get(n,"?")
    if cat_of.get(n)=="QUERY":
        return f"★ {n.replace('_ppk1','')}  →  {REFINE.get(n, place.get(n,{}).get('clade','?'))}"
    if cat_of.get(n)=="OUTGROUP": return f"{name_of.get(n,n)[:28]} [outgroup]"
    return f"{name_of.get(n,n)[:34]} [{c}]"
def lcol(cl):
    return tcol(cl.name) if getattr(cl,"name",None) else "black"

fig=plt.figure(figsize=(15,40)); ax=fig.add_subplot(111)   # tall canvas: 88 tips need vertical room
def supp(c):  # SH-aLRT/UFBoot, on internal branches only (mid-branch → clear of tip labels)
    return c.name if (not c.is_terminal() and c.name and "/" in str(c.name)) else None
Phylo.draw(tree,axes=ax,do_show=False,label_func=lab,label_colors=lcol,
           show_confidence=False,branch_labels=supp)
# style: highlight our sample tips (soothing blue) + shrink branch-support labels
for _t in ax.texts:
    s=_t.get_text().strip()
    if s.startswith("★"):                                   # our sample MAGs
        _t.set_color("#12355B"); _t.set_fontweight("bold"); _t.set_fontsize(6.5)
        _t.set_bbox(dict(boxstyle="round,pad=0.2", facecolor="#CADDF2",
                         edgecolor="#6f9fd8", linewidth=0.7, alpha=0.95))
    elif "/" in s and s.replace("/","").replace(".","").isdigit():   # branch support
        _t.set_fontsize(3.2); _t.set_color("#9a9a9a"); _t.set_va("bottom")
    else:                                                   # reference tip labels
        _t.set_fontsize(4.5)
ax.set_title("Accumulibacter ppk1 ML tree (IQ-TREE TVM+F+I+G4) — placement of 5 CAN MAG ppk1 genes\n"
 "★=our sample MAG (blue highlight); refs blue=Type I, green=Type II, orange=unclassified; rooted on Dechloromonas/Rhodocyclus outgroup",fontsize=11)
ax.set_xlabel("substitutions / site")
for s in ("top","right"): ax.spines[s].set_visible(False)
leg=[Line2D([0],[0],color="#1f77b4",lw=3,label="Type I clades (IA–IE)"),
     Line2D([0],[0],color="#2ca02c",lw=3,label="Type II clades (IIA–IIH)"),
     Line2D([0],[0],color="#ff7f0e",lw=3,label="unclassified DB refs"),
     Line2D([0],[0],marker="s",markerfacecolor="#CADDF2",markeredgecolor="#6f9fd8",
            color="#6f9fd8",lw=0,markersize=12,label="our sample MAG (★, highlighted)")]
ax.legend(handles=leg,loc="upper left",bbox_to_anchor=(1.005,1.0),fontsize=9,frameon=True,
          title="branch support = SH-aLRT / UFBoot")
plt.tight_layout(); plt.savefig(OUTD/"ppk1_reference_tree.png",dpi=150,bbox_inches="tight")
plt.savefig(OUTD/"ppk1_reference_tree.pdf",bbox_inches="tight")
print(f"wrote {OUTD/'ppk1_reference_tree.png'} / .pdf")
