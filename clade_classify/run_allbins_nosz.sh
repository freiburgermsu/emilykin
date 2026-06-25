#!/usr/bin/env bash
# All-bins nosZ pipeline (extraction -> tree -> Clade I/II/III + copy number + server BED).
# Designed to run in a DETACHED tmux session so it survives closing Claude Code / the
# terminal (machine must stay on). Logs (unbuffered) to clade_classify/out/allbins_nosz_run.log.
cd /home/freiburger/Documents/EmilyKin || exit 1
PY="$HOME/Documents/py_venv/bin/python"
LOG="clade_classify/out/allbins_nosz_run.log"
exec > >(tee -a "$LOG") 2>&1     # tee everything to the log AND the tmux pane

echo "=================================================================="
echo "[$(date)] START all-bins nosZ pipeline (tmux session: nosz)"
echo "[$(date)] step 1/2: extract_all_nosz.py"
if ! "$PY" -u clade_classify/extract_all_nosz.py; then
    echo "[$(date)] *** EXTRACT FAILED ***"; exit 1
fi
echo "[$(date)] step 2/2: allbins_tree_classify.py"
if ! "$PY" -u clade_classify/allbins_tree_classify.py; then
    echo "[$(date)] *** CLASSIFY FAILED ***"; exit 1
fi
echo "[$(date)] ALL DONE -> clade_classify/out/{allbins_nosz_clades.tsv,allbins_nosz_per_bin.tsv,"
echo "             allbins_nosz_loci_clade.bed} and clade_classify/server_processing/"
echo "=================================================================="
