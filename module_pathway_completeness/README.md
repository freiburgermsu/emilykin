# Average pathway completeness by co-occurrence module

Per-phase heatmaps of the average metabolic-pathway completeness within each
Louvain module of the per-phase co-occurrence network.

- **Rows** — the 8 pathways summarised in `../ko_pathway_summary.csv` (the pathways
  behind the gene-abundance analysis): Denitrification, DNRA, PolyP, Phosphate,
  PHA, Glycogen, Acetate, Propionate/PHV.
- **Columns** — the Louvain modules defined in that phase
  (`../network/network_module_membership_p_value_FDR_phase{N}.json`); `n` = number
  of MAGs in the module.
- **Cell** — mean, over the organisms in the module, of that pathway's completeness
  (`genes_present / genes_total`). Each module member iterativeID maps 1:1 to a MAG
  via `../mag_iterativeID_old_to_new.json`, and that MAG's per-pathway completeness
  comes from `ko_pathway_summary.csv`.

Files: `module_pathway_completeness_phase{N}.png` (heatmap) and `.csv` (matrix).

**Phase II is absent** — its per-phase co-occurrence network had no FDR-passing
edges, so no Louvain modules are defined. Figures exist for phases I, III, IV, V.

Regenerate: `~/Documents/py_venv/bin/python scripts/render_module_pathway_completeness.py`
