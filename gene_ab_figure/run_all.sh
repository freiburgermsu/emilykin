#!/bin/bash
# Run the full pipeline end-to-end.
# Assumes alignment BAMs already exist in data/ (from 02_build_and_align.sh).
set -euo pipefail

WORK=/scratch1/afreiburger/emilykin/gene_ab_figure
PYBIN=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin/python3

echo "[$(date +%T)] Step 1: Select MAGs and target genes"
cd $WORK && $PYBIN 01_select_mags_and_genes.py

echo "[$(date +%T)] Step 3: Count reads and compute RPKM"
$PYBIN 03_count_and_rpkm.py

echo "[$(date +%T)] Step 4: Plot heatmap"
$PYBIN 04_plot_heatmap.py

echo "[$(date +%T)] Step 5: Bundle offline files"
$PYBIN 05_bundle_offline.py

echo "[$(date +%T)] Done! See gene_abundance_figure.png and offline_bundle/"
