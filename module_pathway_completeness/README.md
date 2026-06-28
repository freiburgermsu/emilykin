# Average pathway completeness by co-occurrence module

Per-phase heatmaps of the average metabolic-pathway completeness within each
Louvain module of the per-phase co-occurrence network.

- **Rows** — 11 rows: the **4 denitrification steps** (see below) followed by the
  7 other pathways summarised in `../ko_pathway_summary.csv` (the pathways behind
  the gene-abundance analysis): DNRA, PolyP, Phosphate, PHA, Glycogen, Acetate,
  Propionate/PHV.
- **Columns** — the Louvain modules defined in that phase
  (`../network/network_module_membership_p_value_FDR_phase{N}.json`); `n` = number
  of MAGs in the module.
- **Cell** — mean, over the organisms in the module, of that row's completeness
  (`genes_present / genes_total`). Each module member iterativeID maps 1:1 to a MAG
  via `../mag_iterativeID_old_to_new.json`, and that MAG's per-pathway completeness
  comes from `ko_pathway_summary.csv`.

## Denitrification disaggregated by step

Instead of a single Denitrification row, the pathway is split into its four
enzymatic steps (reading down: NO₃→NO₂→NO→N₂O→N₂). The four steps partition the
exact same 10 genes that made up the old Denitrification row (5 + 2 + 2 + 1 = 10),
so they are a faithful decomposition of it — the gene-count-weighted mean of the
four step values reproduces the old aggregate per module.

| Row | Step | Genes (gene-count) |
| --- | --- | --- |
| `NO₃→NO₂ (nap/nar)` | nitrate reductase | napA, napB, narG/narZ/nxrA, narH/narY/nxrB, narI/narV (5) |
| `NO₂→NO (nir)` | nitrite reductase | nirK, nirS (2) |
| `NO→N₂O (nor)` | NO reductase | norB, norC (2) |
| `N₂O→N₂ (nos)` | N₂O reductase | nosZ (1) |

The NO₃→NO₂ row combines the periplasmic (**nap**) and membrane-bound (**nar**)
nitrate reductases. Caveat: narG/narH (K00370/K00371) are shared with the
nitrite-oxidoreductase (**nxr**) of nitrifiers, so a high value there can reflect
nitrite oxidisers rather than denitrifiers.

Files: `module_pathway_completeness_phase{N}.png` (heatmap) and `.csv` (matrix).

**Phase II is absent** — its per-phase co-occurrence network had no FDR-passing
edges, so no Louvain modules are defined. Figures exist for phases I, III, IV, V.

Regenerate: `~/Documents/py_venv/bin/python scripts/render_module_pathway_completeness.py`
