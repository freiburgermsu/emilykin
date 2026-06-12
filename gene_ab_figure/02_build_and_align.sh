#!/bin/bash
# Step 2: Build combined MAG reference FASTA (unique contig names) and align long reads

set -euo pipefail

WORK=/scratch1/afreiburger/emilykin/gene_ab_figure
DREP=/scratch1/afreiburger/emilykin/processed/mag/drep/dereplicated_genomes
READS=/scratch1/afreiburger/emilykin/raw/meta/longreads
DATA=$WORK/data

MINIMAP2=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/d501a828539edb92b78db288661dd07f_/bin/minimap2
SAMTOOLS=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/88fdb48d4d745c55ec2cd90b407de422_/bin/samtools

THREADS=16
REF=$DATA/combined_ref.fa
BAM=$DATA/aligned_sorted.bam

# ── 1. Build combined reference FASTA ────────────────────────────────────
echo "Building combined reference FASTA with unique contig names..."
rm -f $REF
while IFS= read -r mag; do
    fa=$DREP/${mag}.fa
    if [[ ! -f $fa ]]; then
        echo "  WARNING: $fa not found, skipping"
        continue
    fi
    safe_mag="${mag//./_}"
    # Prefix each contig header with MAG name
    awk -v prefix="${safe_mag}::" '
        /^>/ { print ">" prefix substr($0, 2); next }
        { print }
    ' $fa >> $REF
done < $DATA/selected_mag_list.txt

echo "  Combined reference: $(grep -c "^>" $REF) contigs"
echo "  Size: $(du -sh $REF | cut -f1)"

# ── 2. Index reference ─────────────────────────────────────────────────
echo "Indexing reference with samtools faidx..."
$SAMTOOLS faidx $REF

# ── 3. Align all long reads (all samples combined) ─────────────────────
echo "Aligning long reads from all samples..."
ALL_READS=$(ls $READS/*.fastq.gz | tr '\n' ' ')
echo "  Read files: $(ls $READS/*.fastq.gz | wc -l)"

$MINIMAP2 \
    -ax map-ont \
    -t $THREADS \
    --secondary=no \
    $REF \
    $ALL_READS \
| $SAMTOOLS sort -@ $THREADS -o $BAM -

$SAMTOOLS index $BAM

echo "  Aligned BAM: $BAM"
echo "  Total reads mapped: $($SAMTOOLS flagstat $BAM | grep 'mapped (' | head -1)"

# ── 4. Also align per-sample for per-sample RPKM if needed ──────────────
for sample in CAN_1 CAN_2 CAN_3 CAN_4 CAN_5; do
    reads=$READS/${sample}_nanopore.fastq.gz
    bam_out=$DATA/aligned_${sample}_sorted.bam
    if [[ ! -f $reads ]]; then
        echo "  WARNING: $reads not found"
        continue
    fi
    echo "  Aligning $sample..."
    $MINIMAP2 -ax map-ont -t $THREADS --secondary=no $REF $reads \
    | $SAMTOOLS sort -@ $THREADS -o $bam_out -
    $SAMTOOLS index $bam_out
    echo "  $sample mapped: $($SAMTOOLS flagstat $bam_out | grep 'mapped (' | head -1)"
done

echo "Done. Next: run 03_count_and_rpkm.py"
