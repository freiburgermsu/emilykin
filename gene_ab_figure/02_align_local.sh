#!/bin/bash
# Step 2 (local): align each long-read sample to the 276-genome community
# reference, streaming minimap2 PAF directly into count_stream.py. No samtools,
# no on-disk PAF/BAM. Produces data/counts_<sample>.json per sample.
set -uo pipefail

WORK=/home/freiburger/Documents/EmilyKin/gene_ab_figure
DATA=$WORK/data
READS=/home/freiburger/Documents/EmilyKin/meta/longreads
MM2=/home/freiburger/Documents/EmilyKin/meta/mag_gene_ab/tools/minimap2-2.28_x64-linux/minimap2
PY=/home/freiburger/Documents/py_venv/bin/python
REF=$DATA/combined_ref.mmi
BED=$DATA/target_genes.bed
THREADS=48

mkdir -p $DATA/logs
echo "[$(date +%T)] START alignment of 5 long-read samples (ref=276 MAGs)"

for sample in CAN_1 CAN_2 CAN_3 CAN_4 CAN_5; do
    reads=$READS/${sample}_nanopore.fastq.gz
    out=$DATA/counts_${sample}.json
    log=$DATA/logs/minimap2_${sample}.log
    if [[ ! -f $reads ]]; then echo "  MISSING $reads — skipping"; continue; fi
    if [[ -f $out ]]; then echo "[$(date +%T)] $sample already counted — skipping"; continue; fi
    echo "[$(date +%T)] aligning $sample ($(du -h $reads | cut -f1)) ..."
    $MM2 -x map-ont -t $THREADS --secondary=no -K 2g "$REF" "$reads" 2> "$log" \
        | $PY $WORK/count_stream.py --bed "$BED" --sample "$sample" --out "$out"
    rc=${PIPESTATUS[0]}
    echo "[$(date +%T)] $sample done (minimap2 rc=$rc) -> $out"
done

echo "[$(date +%T)] ALL SAMPLES DONE"
touch $DATA/.alignment_complete
