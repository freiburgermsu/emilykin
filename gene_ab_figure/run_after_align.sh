#!/bin/bash
# Waits for the tmux alignment to finish, then aggregates RPKM and renders the figure.
set -uo pipefail
cd /home/freiburger/Documents/EmilyKin/gene_ab_figure
PY=/home/freiburger/Documents/py_venv/bin/python
while [ ! -f data/.alignment_complete ]; do sleep 120; done
echo "[$(date +%T)] alignment complete — aggregating RPKM"
$PY 03_aggregate_rpkm_local.py
echo "[$(date +%T)] rendering figure"
$PY 04_plot_heatmap.py
echo "[$(date +%T)] FIGURE DONE"
