#!/usr/bin/env bash
# Build a 276-MAG dereplicated combined reference (::-naming, matching the
# allbins_nosz BED + the focal-BAM @SQ style) and map the 5 Nanopore samples to
# it, so compute_nosz_rpkm.py can quantify ALL nosZ-bearing bins (not just the
# 21 focal MAGs whose BAMs already existed).
set -euo pipefail

MM2=/scratch/afreiburger/metag-hybrid/env/miniforge/envs/assembly/bin/minimap2
SAM=/scratch/afreiburger/metag-hybrid/env/miniforge/envs/assembly/bin/samtools
MAGDIR=/scratch1/afreiburger/emilykin/processed/mag/drep/dereplicated_genomes
READS=/scratch1/afreiburger/emilykin/processed/reads/nanopore
WORK=/scratch/afreiburger/nosz_allbins          # local NVMe
TMP=$WORK/tmp
REF=$WORK/derep276_ref.fa
MMI=$WORK/derep276_ref.mmi
LOG=$WORK/build.log
SAMPLES=(CAN_1 CAN_2 CAN_3 CAN_4 CAN_5)
mkdir -p "$WORK" "$TMP"

echo "[$(date -u +%FT%TZ)] START build_and_map_276ref" | tee -a "$LOG"

# 1. Build prefixed reference: >contig -> >MAG_underscores::contig
if [[ ! -f "$REF" ]]; then
  echo "[$(date -u +%FT%TZ)] building reference from $(ls "$MAGDIR"/*.fa | wc -l) MAGs" | tee -a "$LOG"
  for fa in "$MAGDIR"/*.fa; do
    mag=$(basename "$fa" .fa); safe=${mag//./_}
    awk -v m="$safe" '/^>/{split($1,a,">"); print ">" m "::" a[2]; next}{print}' "$fa"
  done > "$REF"
fi
NCTG=$(grep -c '^>' "$REF"); NDUP=$(grep '^>' "$REF" | sort | uniq -d | wc -l)
echo "[$(date -u +%FT%TZ)] reference: $NCTG contigs, $NDUP duplicate names" | tee -a "$LOG"

# 2. minimap2 index (map-ont preset)
if [[ ! -f "$MMI" ]]; then
  echo "[$(date -u +%FT%TZ)] indexing -> $MMI" | tee -a "$LOG"
  "$MM2" -x map-ont -d "$MMI" "$REF" 2>>"$LOG"
fi

# 3. map each Nanopore sample (2 concurrent), primary-only sorted BAM
map_one() {
  local s="$1"
  local bam="$WORK/aligned_${s}_sorted.bam"
  if [[ -f "${bam}.bai" ]]; then echo "[$(date -u +%FT%TZ)] $s already done, skip" >>"$LOG"; return; fi
  echo "[$(date -u +%FT%TZ)] mapping $s" >>"$LOG"
  "$MM2" -ax map-ont -t 48 --secondary=no "$MMI" "$READS/${s}.filt.fastq.gz" 2>>"$LOG" \
    | "$SAM" sort -@ 6 -m 3G -T "$TMP/${s}_sort" -o "$bam" -
  "$SAM" index "$bam"
  local n; n=$("$SAM" idxstats "$bam" | awk '{t+=$3} END{print t}')
  echo "[$(date -u +%FT%TZ)] $s DONE: ${n} primary-mapped reads" >>"$LOG"
}
export -f map_one; export MM2 SAM WORK READS MMI TMP LOG

printf '%s\n' "${SAMPLES[@]}" | xargs -P 2 -I{} bash -c 'map_one "$@"' _ {}

echo "[$(date -u +%FT%TZ)] ALL MAPPING DONE" | tee -a "$LOG"
ls -la "$WORK"/aligned_*_sorted.bam | tee -a "$LOG"
