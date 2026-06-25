#!/usr/bin/env python3
"""
Part D — definitive tree-based clade I/II placement (resolves CLADE_I_II_DISCREPANCY.md).

Reference+query ML tree: out/nosz_tree.nwk (FastTree WAG+GAMMA, SH-like local
supports) built from 330 Chee+Orellana C-NosZ refs (clade I+II, genus-named),
259 He 269NosZ refs, and the 13 K00376 query proteins.

Independent method (does NOT use the `1577*` file-name clade labels or diamond %id):
  1. Label reference tips by clade with ONLY textbook-unambiguous marker genera.
     The disputed Rhodocyclaceae/Comamonadaceae (Dechloromonas, Dechlorosoma,
     Accumulibacter, Azonexus, Alicycliphilus, Acidovorax, Magnetospirillum) are
     NOT used as anchors, so the tree places them independently.
  2. Clade II clan = MRCA(clade-II anchors). It is a single-edge bipartition
     (rooting-independent) that must contain every clade-II anchor and zero
     clade-I anchors. SH support of that edge is reported.
  3. A tip is Clade II iff it lies inside that clan; otherwise (for a C-NosZ tip)
     Clade I — corroborated by its sister reference clan (climb to the smallest
     clade holding a labelled reference) and that clade's SH support.
  4. Re-examine NosZREF-1577D: report the clan each 1577D reference taxon falls in.
  5. Reconcile vs the HMM (figure) and diamond calls read from the TSVs.

Writes out/partD_tree_clades.tsv, out/partD_reference_audit.tsv.
"""
from pathlib import Path
import csv
from collections import Counter
from Bio import Phylo

HERE = Path(__file__).parent
OUT  = HERE / "out"
TREE = OUT / "nosz_tree.nwk"

# textbook-unambiguous anchors (genus substrings); disputed taxa excluded on purpose
I_ANCH = ["Pseudomonas","Paracoccus","Bradyrhizobium","Sinorhizobium","Rhodobacter",
          "Rhodopseudomonas","Shewanella","Marinobacter","Brucella","Azoarcus",
          "Aromatoleum","Thauera","Achromobacter","Ralstonia","Cupriavidus",
          "Neisseria","Kingella","Roseobacter","Halomonas"]
II_ANCH = ["Wolinella","Anaeromyxobacter","Gemmatimonas","Chryseobacterium","Prevotella",
           "Salinibacter","Rhodothermus","Flavobacterium","Cytophaga","Dyadobacter",
           "Haliscomenobacter","Gramella","Kordia","Persicobacter","Chlorobaculum",
           "Ignavibacterium","Bacteroides","Sphingobacterium"]
DISPUTED = ["Dechloromonas","Dechlorosoma","Accumulibacter","Azonexus","Alicycliphilus",
            "Acidovorax","Magnetospirillum"]

QUERY_MAG = {
    "CAN_1_bin_77_DLHCAE_04715":"CAN_1_bin.77",  "CAN_1_bin_98_KLJHLB_03109":"CAN_1_bin.98",
    "CAN_2_bin_6_HMFKDB_00723":"CAN_2_bin.6",    "CAN_3_bin_203_NGIEFB_03466":"CAN_3_bin.203",
    "CAN_3_bin_221_ECOIAL_03204":"CAN_3_bin.221","CAN_5_bin_112_DCCNNJ_01028":"CAN_5_bin.112",
    "CAN_5_bin_147_CIIJDE_02562":"CAN_5_bin.147","CAN_5_bin_40_GNOLDL_01230":"CAN_5_bin.40",
    "CAN_5_bin_70_MDDBPK_01672":"CAN_5_bin.70",  "coasm_bin_185_FDDLBN_03425":"coasm_bin.185",
    "coasm_bin_260_FLKOMO_04520":"coasm_bin.260","coasm_bin_481_EGDGHF_00332":"coasm_bin.481",
    "coasm_bin_55_FCMEHI_03805":"coasm_bin.55",
}
MAG2GID = {  # for joining to the TSVs (gene_id) and reporting
    "CAN_1_bin.77":"CAN_1_bin.77|DLHCAE_04715","CAN_1_bin.98":"CAN_1_bin.98|KLJHLB_03109",
    "CAN_2_bin.6":"CAN_2_bin.6|HMFKDB_00723","CAN_3_bin.203":"CAN_3_bin.203|NGIEFB_03466",
    "CAN_3_bin.221":"CAN_3_bin.221|ECOIAL_03204","CAN_5_bin.112":"CAN_5_bin.112|DCCNNJ_01028",
    "CAN_5_bin.147":"CAN_5_bin.147|CIIJDE_02562","CAN_5_bin.40":"CAN_5_bin.40|GNOLDL_01230",
    "CAN_5_bin.70":"CAN_5_bin.70|MDDBPK_01672","coasm_bin.185":"coasm_bin.185|FDDLBN_03425",
    "coasm_bin.260":"coasm_bin.260|FLKOMO_04520","coasm_bin.481":"coasm_bin.481|EGDGHF_00332",
    "coasm_bin.55":"coasm_bin.55|FCMEHI_03805",
}

def is_ref(n):  return n.startswith("REF__")
def is_he(n):   return n.startswith("HE__")
def is_query(n): return not (is_ref(n) or is_he(n))
def matches(n, keys): return any(k.lower() in n.lower() for k in keys)
def genus(n): return n.replace("REF__","").split("_")[0]

# ── prior calls for comparison ──────────────────────────────────────────────
def read_clades(path, col):
    d = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[r["mag"]] = r[col]
    return d
hmm     = read_clades("../gene_ab_figure/data/nosz_clades.tsv", "clade")   # figure (HMM)
diamond = read_clades(OUT / "nosz_clades_updated.tsv", "clade")            # diamond
subclade = read_clades("../gene_ab_figure/data/nosz_clades.tsv", "cladeII_best_subclade")

# ── tree + robust Clade II clan (midpoint root → single-edge bipartition) ────
tree = Phylo.read(TREE, "newick")
tree.root_at_midpoint()
tips = [l.name for l in tree.get_terminals()]
I_tips  = [t for t in tips if is_ref(t) and matches(t, I_ANCH)]
II_tips = [t for t in tips if is_ref(t) and matches(t, II_ANCH)]
mII = tree.common_ancestor([l for l in tree.get_terminals() if l.name in II_tips])
sII = {l.name for l in mII.get_terminals()}                 # the Clade II clan
n_I_in_II  = sum(t in set(I_tips) for t in sII)
n_II_in_II = sum(t in set(II_tips) for t in sII)
print(f"Clade II clan: {len(sII)} tips | clade-II anchors inside = {n_II_in_II}/{len(II_tips)}"
      f" | clade-I anchors inside = {n_I_in_II}/{len(I_tips)} | SH support = {mII.confidence}")
assert n_I_in_II == 0 and n_II_in_II == len(II_tips), "Clade II clan not clean!"

def call_tip(name):           # robust: inside the clade II clan ⇒ II, else I (for C-NosZ)
    return "II" if name in sII else "I"

name2leaf = {l.name: l for l in tree.get_terminals()}
def sister(qname):
    """Smallest clade containing the query and ≥1 labelled reference."""
    path = tree.get_path(name2leaf[qname])
    for node in reversed([tree.root] + path[:-1]):
        refs = [t.name for t in node.get_terminals() if not is_query(t.name)]
        if refs:
            calls = Counter(call_tip(m) for m in refs)
            gen = [genus(m) for m in refs if is_ref(m)]
            return calls, node.confidence, [g for g, _ in Counter(gen).most_common(5)]
    return Counter(), None, []

# ── per-query placement ─────────────────────────────────────────────────────
print("\n=== TREE-BASED QUERY PLACEMENT vs prior calls ===")
hdr = f"{'MAG':14s} {'HMM':4s} {'DIA':4s} {'TREE':5s} {'supp':6s} sister(clade:n) / nearest clade-I/II ref genera"
print(hdr); print("-"*len(hdr))
rows = []
for q, mag in sorted(QUERY_MAG.items(), key=lambda kv: kv[1]):
    tree_clade = call_tip(q)
    calls, supp, gen = sister(q)
    rows.append({
        "mag": mag, "gene_id": MAG2GID[mag], "tree_clade": tree_clade,
        "placement_support": "" if supp is None else f"{supp:.3f}",
        "sister_composition": dict(calls), "nearest_ref_genera": ";".join(gen),
        "hmm_clade": hmm.get(mag,"?"), "diamond_clade": diamond.get(mag,"?"),
        "cladeII_subclade": subclade.get(mag,""),
    })
    print(f"{mag:14s} {hmm.get(mag,'?'):4s} {diamond.get(mag,'?'):4s} {tree_clade:5s} "
          f"{('' if supp is None else f'{supp:.3f}'):6s} {dict(calls)}  {gen}")

# ── step 2: NosZREF-1577D reference taxa — which clan? ───────────────────────
print("\n=== STEP 2 — disputed reference taxa (NosZREF-1577D members) clan membership ===")
audit = []
for t in tips:
    if is_ref(t) and matches(t, DISPUTED):
        audit.append({"genus": genus(t), "clade": call_tip(t), "tip": t})
by_gen = {}
for a in audit:
    by_gen.setdefault(a["genus"], Counter())[a["clade"]] += 1
for g, c in sorted(by_gen.items()):
    print(f"  {g:18s} {dict(c)}")

# ── summary / reconciliation ────────────────────────────────────────────────
nI = sum(1 for r in rows if r["tree_clade"] == "I")
print(f"\nTREE RESULT: Clade I = {nI}   Clade II = {len(rows)-nI}")
agree_hmm = sum(1 for r in rows if r["tree_clade"] == r["hmm_clade"])
agree_dia = sum(1 for r in rows if r["tree_clade"] == r["diamond_clade"])
print(f"  agree with HMM(figure): {agree_hmm}/13   agree with diamond: {agree_dia}/13")
print("  TREE vs DIAMOND differences:",
      [r["mag"] for r in rows if r["tree_clade"] != r["diamond_clade"]])
print("  TREE vs HMM differences:",
      [r["mag"] for r in rows if r["tree_clade"] != r["hmm_clade"]])

# ── write ───────────────────────────────────────────────────────────────────
with open(OUT / "partD_tree_clades.tsv", "w", newline="") as f:
    cols = ["mag","gene_id","tree_clade","placement_support","sister_composition",
            "nearest_ref_genera","hmm_clade","diamond_clade","cladeII_subclade"]
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader()
    for r in rows:
        r = dict(r); r["sister_composition"] = str(r["sister_composition"]); w.writerow(r)
with open(OUT / "partD_reference_audit.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["genus","clade","tip"], delimiter="\t")
    w.writeheader(); w.writerows(audit)
print(f"\nwrote {OUT}/partD_tree_clades.tsv  and  partD_reference_audit.tsv")
