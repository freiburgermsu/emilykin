#!/bin/bash
# Align the 5 Nanopore samples to the 276-MAG community reference and count reads
# over the KO-gene BED (all 42 KOs). -ax map-ont (base-level) so reads land on the
# correct homologous gene copy; SAM streamed into count_stream_sam.py. Same method
# as the (validated) gene_ab_figure alignment, extended to the full 276-MAG community.
set -uo pipefail

REPO=/home/freiburger/Documents/EmilyKin
WORK=$REPO/ko_gene_abundance
READS=$REPO/meta/longreads
MM2=$REPO/meta/mag_gene_ab/tools/minimap2-2.28_x64-linux/minimap2
PY=/home/freiburger/Documents/py_venv/bin/python
REF=$REPO/gene_ab_figure/data/combined_ref.mmi          # 276-MAG community index
BED=$WORK/ko_genes.bed
COUNTER=$REPO/gene_ab_figure/count_stream_sam.py
THREADS=48

mkdir -p $WORK/logs
rm -f $WORK/.ko_align_complete
echo "[$(date +%T)] START KO alignment (276-MAG ref; $(wc -l <$BED) BED rows)"
for s in CAN_1 CAN_2 CAN_3 CAN_4 CAN_5; do
    reads=$READS/${s}_nanopore.fastq.gz
    out=$WORK/counts_${s}.json
    log=$WORK/logs/mm2_${s}.log
    [[ -f $reads ]] || { echo "  MISSING $reads"; continue; }
    echo "[$(date +%T)] aligning $s ($(du -h $reads | cut -f1)) ..."
    $MM2 -ax map-ont -t $THREADS --secondary=no -K 2g "$REF" "$reads" 2>"$log" \
        | $PY $COUNTER --bed "$BED" --sample "$s" --out "$out"
    echo "[$(date +%T)] $s done (rc=${PIPESTATUS[0]}) -> $out"
done
echo "[$(date +%T)] ALL DONE"
touch $WORK/.ko_align_complete
