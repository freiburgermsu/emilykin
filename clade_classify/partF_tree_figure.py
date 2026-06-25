#!/usr/bin/env python3
"""
Part F — emit the taxonomical reference tree as inspectable artifacts:
  out/nosz_tree_annotated.nwk : full 602-tip ML tree, every tip relabelled with its
                                clade ([I]/[II]/[III]) and queries marked [Q>I]/[Q>II].
  out/nosz_reference_tree.png : pruned representative tree (clade-I/II/III anchors,
                                the disputed Rhodocyclaceae, + all 13 queries) showing
                                the placements that resolve the discrepancy.
"""
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

HERE = Path(__file__).parent; OUT = HERE / "out"
I_ANCH = ["Pseudomonas","Paracoccus","Bradyrhizobium","Sinorhizobium","Rhodobacter",
          "Rhodopseudomonas","Shewanella","Marinobacter","Brucella","Azoarcus",
          "Aromatoleum","Thauera","Achromobacter","Ralstonia","Cupriavidus",
          "Neisseria","Kingella","Roseobacter","Halomonas"]
II_ANCH = ["Wolinella","Anaeromyxobacter","Gemmatimonas","Chryseobacterium","Prevotella",
           "Salinibacter","Rhodothermus","Flavobacterium","Cytophaga","Dyadobacter",
           "Haliscomenobacter","Gramella","Kordia","Persicobacter","Chlorobaculum",
           "Ignavibacterium","Bacteroides","Sphingobacterium"]
def is_ref(n):  return n.startswith("REF__")
def is_he(n):   return n.startswith("HE__")
def is_query(n): return not (is_ref(n) or is_he(n))
def matches(n, keys): return any(k.lower() in n.lower() for k in keys)
def genus(n): return n.replace("REF__","").replace("HE__","").split("_")[0]

tree = Phylo.read(OUT / "nosz_tree.nwk", "newick")
tree.root_at_midpoint()
tips = [l.name for l in tree.get_terminals()]
I_tips  = [t for t in tips if is_ref(t) and matches(t, I_ANCH)]
II_tips = [t for t in tips if is_ref(t) and matches(t, II_ANCH)]
sII  = {l.name for l in tree.common_ancestor([l for l in tree.get_terminals() if l.name in II_tips]).get_terminals()}
s545 = {l.name for l in tree.common_ancestor([l for l in tree.get_terminals() if l.name in I_tips]).get_terminals()}
def clade_of(n):                      # I / II / III(out)
    if n in sII:  return "II"
    if n in s545: return "I"
    return "III"

# ── 1. annotated full Newick ────────────────────────────────────────────────
t_full = Phylo.read(OUT / "nosz_tree.nwk", "newick")
for l in t_full.get_terminals():
    c = clade_of(l.name)
    tag = (f"Q>{c}" if is_query(l.name) else c)
    base = l.name if is_query(l.name) else f"{genus(l.name)}"
    l.name = f"[{tag}]{base}"
Phylo.write(t_full, OUT / "nosz_tree_annotated.nwk", "newick")
print(f"wrote {OUT/'nosz_tree_annotated.nwk'}")

# ── 2. pruned representative tree figure ────────────────────────────────────
REPS = {  # representative genera to keep (incl. the disputed Rhodocyclaceae)
 "I":  ["Pseudomonas","Paracoccus","Bradyrhizobium","Shewanella","Azoarcus","Thauera",
        "Ralstonia","Cupriavidus","Rubrivivax","Burkholderia","Rhodoferax","Aromatoleum",
        "Acidovorax","Alicycliphilus"],
 "II": ["Wolinella","Anaeromyxobacter","Gemmatimonas","Chryseobacterium","Flavisolibacter",
        "Dyadobacter","Caldilinea","Rhodothermus","Sphaerobacter","Thermomicrobium",
        "Dechloromonas","Dechlorosoma","Accumulibacter","Magnetospirillum"],
}
keep = set(t for t in tips if is_query(t))
for c, gens in REPS.items():
    for g in gens:
        got = [t for t in tips if is_ref(t) and genus(t).lower() == g.lower() and clade_of(t) == c]
        keep.update(got[:2])
keep.update([t for t in tips if clade_of(t) == "III"][:3])   # a few outgroup tips

for l in list(tree.get_terminals()):
    if l.name not in keep:
        tree.prune(l)
print(f"pruned to {len(tree.get_terminals())} tips for the figure")

COL = {"I":"#1f77b4", "II":"#2ca02c", "III":"#7f7f7f"}
def label(cl):
    if not cl.name: return None
    c = clade_of(cl.name)
    if is_query(cl.name):
        return f"★ {cl.name}  [Clade {c}]"
    return f"{genus(cl.name)} [{c}]"
def lblcol(cl):
    n = getattr(cl, "name", None)
    if not n: return "black"
    return "#12355B" if is_query(n) else COL[clade_of(n)]

def _supp(c):  # FastTree SH-like local support (0-1) on internal branches only
    return f"{c.confidence:.2f}" if (not c.is_terminal() and c.confidence is not None) else None
fig = plt.figure(figsize=(13, 16)); ax = fig.add_subplot(111)
Phylo.draw(tree, axes=ax, do_show=False, label_func=label,
           label_colors=lblcol, show_confidence=False, branch_labels=_supp)
# style: highlight our sample tips (soothing blue) + shrink branch-support labels
for _t in ax.texts:
    s = _t.get_text().strip()
    if s.startswith("★"):
        _t.set_color("#12355B"); _t.set_fontweight("bold")
        _t.set_fontsize(_t.get_fontsize() + 1)
        _t.set_bbox(dict(boxstyle="round,pad=0.3", facecolor="#CADDF2",
                         edgecolor="#6f9fd8", linewidth=0.8, alpha=0.95))
    elif s.replace(".", "").isdigit() and len(s) <= 4:
        _t.set_fontsize(4.2); _t.set_color("#8a8a8a"); _t.set_va("bottom")
ax.set_title("nosZ reference tree (Chee+Orellana C-NosZ + He refs) — tree placement of the 13 CAN K00376 genes\n"
             "★ = our sample MAG (blue highlight); references coloured by clade   blue=I  green=II  grey=III/outgroup",
             fontsize=11)
ax.set_xlabel("substitutions / site");
for sp in ("top","right"): ax.spines[sp].set_visible(False)
leg = [Line2D([0],[0],color=COL["I"],lw=3,label="Clade I (typical, Pseudomonadota)"),
       Line2D([0],[0],color=COL["II"],lw=3,label="Clade II (atypical; incl. Dechloromonas/Accumulibacter)"),
       Line2D([0],[0],color=COL["III"],lw=3,label="Clade III / outgroup (He refs)"),
       Line2D([0],[0],marker="s",markerfacecolor="#CADDF2",markeredgecolor="#6f9fd8",
              color="#6f9fd8",lw=0,markersize=12,label="our sample MAG (★, highlighted)")]
ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=9,
          frameon=True, title="branch support = FastTree SH-like (0–1)")
plt.tight_layout()
plt.savefig(OUT / "nosz_reference_tree.png", dpi=150, bbox_inches="tight")
plt.savefig(OUT / "nosz_reference_tree.pdf", bbox_inches="tight")
print(f"wrote {OUT/'nosz_reference_tree.png'} / .pdf")
