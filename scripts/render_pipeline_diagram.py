#!/usr/bin/env python3
"""Render a multi-page PDF documenting the full emilykin metagenomic pipeline.

Pages
-----
1. Main assembly → binning → MAG QC → annotation pipeline
2. 16S amplicon → MAG abundance bridge  +  direct metagenomic abundance
3. Bakta vs alternative annotation tools (justified comparison table)

Run:
    python render_pipeline_diagram.py  (produces pipeline_diagram.pdf)
"""

import textwrap
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

OUT = "pipeline_diagram.pdf"

# ── colour palette ──────────────────────────────────────────────────────────
C = {
    "qc":     "#4a90d9",   # blue
    "asm":    "#5cb85c",   # green
    "bin":    "#f0ad4e",   # orange
    "mag":    "#9b59b6",   # purple
    "ann":    "#e74c3c",   # red
    "abund":  "#1abc9c",   # teal
    "arrow":  "#555555",
    "bg":     "#f8f9fa",
    "header": "#2c3e50",
    "note":   "#7f8c8d",
    "border": "#ecf0f1",
}

FONT = "DejaVu Sans"

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def box(ax, x, y, w, h, color, label, sublabel=None, fontsize=9, radius=0.015, alpha=0.92):
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        facecolor=color, edgecolor="white", linewidth=1.5, alpha=alpha,
        zorder=3,
    )
    ax.add_patch(patch)
    ty = y + (h * 0.12 if sublabel else 0)
    ax.text(x, ty, label, ha="center", va="center", fontsize=fontsize,
            fontfamily=FONT, fontweight="bold", color="white", zorder=4,
            wrap=True)
    if sublabel:
        ax.text(x, y - h * 0.25, sublabel, ha="center", va="center",
                fontsize=fontsize - 1.5, fontfamily=FONT, color="white",
                alpha=0.88, zorder=4, style="italic",
                multialignment="center")


def arrow(ax, x0, y0, x1, y1, label=None, color=C["arrow"], lw=1.5,
          connectionstyle="arc3,rad=0.0"):
    ax.annotate(
        "", xy=(x1, y1), xytext=(x0, y0),
        arrowprops=dict(
            arrowstyle="-|>", color=color, lw=lw,
            connectionstyle=connectionstyle,
        ),
        zorder=2,
    )
    if label:
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        ax.text(mx + 0.01, my, label, ha="left", va="center",
                fontsize=7, fontfamily=FONT, color=C["note"], zorder=5)


def section_header(ax, x, y, text, color=C["header"], fontsize=11):
    ax.text(x, y, text, ha="left", va="center",
            fontsize=fontsize, fontfamily=FONT, fontweight="bold",
            color=color, zorder=5)
    ax.plot([x, x + 0.88], [y - 0.012, y - 0.012],
            color=color, lw=1, alpha=0.4, zorder=4)


def param_note(ax, x, y, lines, fontsize=6.5, color=C["note"]):
    txt = "\n".join(lines)
    ax.text(x, y, txt, ha="center", va="top",
            fontsize=fontsize, fontfamily=FONT, color=color,
            multialignment="center", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=C["border"], alpha=0.85))


def page_background(fig):
    fig.patch.set_facecolor(C["bg"])


# ─────────────────────────────────────────────────────────────────────────────
# PAGE 1 – Full assembly / binning / annotation pipeline
# ─────────────────────────────────────────────────────────────────────────────

def make_page1(pdf):
    fig, ax = plt.subplots(figsize=(17, 22))
    page_background(fig)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_facecolor(C["bg"])

    # ── Title ──
    ax.text(0.5, 0.975, "EmilyKin Hybrid Metagenomics Pipeline",
            ha="center", va="top", fontsize=16, fontfamily=FONT,
            fontweight="bold", color=C["header"])
    ax.text(0.5, 0.960, "5 samples (CAN 1–5)  ·  Nanopore ONT + Illumina NovaSeq  ·  Snakemake workflow",
            ha="center", va="top", fontsize=9, fontfamily=FONT, color=C["note"])

    # ── Column x centres ──
    xL, xC, xR = 0.18, 0.50, 0.82   # per-sample hybrid | co-assembly | per-sample short

    # ── Column labels ──
    for xc, lbl in [(xL, "Per-Sample Hybrid\n(×5 samples)"),
                     (xC, "Hybrid Co-Assembly\n(all 5 concatenated)"),
                     (xR, "Per-Sample Short\n(×5, metaSPAdes)")]:
        ax.text(xc, 0.940, lbl, ha="center", va="top", fontsize=8.5,
                fontfamily=FONT, fontweight="bold", color=C["asm"],
                multialignment="center")

    bw, bh = 0.28, 0.045   # box width, height

    # ── ROW: Raw reads ──
    y = 0.885
    box(ax, xL, y, bw, bh, C["qc"],
        "Raw ONT reads",
        "filtlong  min_len=1 kb  keep=95 %",
        fontsize=8)
    box(ax, xC, y, bw, bh, C["qc"],
        "Pooled ONT reads (CAN 1–5 cat)",
        "filtlong  min_len=1 kb  keep=95 %",
        fontsize=8)
    box(ax, xR, y, bw, bh, C["qc"],
        "Raw Illumina PE reads",
        "fastp  Q20 tail trim  len≥50  PE corr.",
        fontsize=8)

    # Illumina QC also feeds left & centre (after fastp)
    for xc in (xL, xC):
        ax.annotate("", xy=(xc - 0.02, y + bh/2), xytext=(xR - bw/2 - 0.02, y),
                    arrowprops=dict(arrowstyle="-|>", color=C["qc"], lw=1,
                                   connectionstyle="arc3,rad=-0.2"), zorder=2)
    ax.text(0.375, y + 0.025, "fastp clean reads\n→ per-sample + co-assembly",
            ha="center", va="bottom", fontsize=6.5, color=C["note"],
            fontfamily=FONT, style="italic")

    # ── STEP 1: Flye ──
    y1 = 0.810
    for xc, lbl in [(xL, "metaFlye  (per-sample)"),
                    (xC, "metaFlye  (co-assembly)")]:
        arrow(ax, xc, y - bh/2, xc, y1 + bh/2)
        box(ax, xc, y1, bw, bh, C["asm"], lbl,
            "--meta  --nano-hq  threads=64",
            fontsize=8)
    param_note(ax, xL - 0.20, y1 + bh/2 - 0.01,
               ["Assembly", "strategy:","De Bruijn graph",
                "+ repeat graph", "for metagenomes"],
               fontsize=6)

    # ── STEP 2: Medaka ──
    y2 = 0.735
    for xc in (xL, xC):
        arrow(ax, xc, y1 - bh/2, xc, y2 + bh/2)
        box(ax, xc, y2, bw, bh, C["asm"],
            "Medaka  (ONT neural-net polish)",
            "model: r1041_e82_400bps_sup_v5.0.0\nbatch=200  2×H100 GPU parallel",
            fontsize=7.5)
    param_note(ax, xL - 0.20, y2 + bh/2 - 0.01,
               ["Contigs split","across GPU 0","and GPU 1","(mini_align)"],
               fontsize=6)

    # ── STEP 3: Polypolish ──
    y3 = 0.660
    for xc in (xL, xC):
        arrow(ax, xc, y2 - bh/2, xc, y3 + bh/2)
        box(ax, xc, y3, bw, bh, C["asm"],
            "Polypolish  (Illumina short-read polish)",
            "bwa mem -a  →  polypolish filter  →  polish",
            fontsize=8)
    param_note(ax, xL - 0.20, y3 + bh/2 - 0.01,
               ["Corrects SNPs","& indels using","per-read","alignments"],
               fontsize=6)

    # ── STEP 4: Pilon ──
    y4 = 0.585
    for xc in (xL, xC):
        arrow(ax, xc, y3 - bh/2, xc, y4 + bh/2)
        box(ax, xc, y4, bw, bh, C["asm"],
            "Pilon  (final Illumina polish)",
            "bwa-mem2  →  pilon --fix snps,indels\n-Xmx256g  threads=32",
            fontsize=8)

    # ── metaSPAdes (right column) ──
    arrow(ax, xR, y - bh/2, xR, y4 + bh/2)
    box(ax, xR, y4, bw, bh, C["asm"],
        "metaSPAdes  (per-sample short-read)",
        "--meta  -k 21,33,55,77,99,127\nmem=600 GB  threads=64",
        fontsize=8)

    # ── "final.fasta" convergence ──
    y4b = 0.533
    for xc in (xL, xC, xR):
        arrow(ax, xc, y4 - bh/2, xc, y4b + 0.008)
    ax.text(0.50, y4b, "▼  Final assemblies per source  ▼",
            ha="center", va="center", fontsize=8, fontfamily=FONT,
            color=C["note"], style="italic")

    # ─── BINNING section ───────────────────────────────────────────────────
    section_header(ax, 0.03, 0.508, "BINNING  (MetaBAT2 — hybrid assemblies + co-assembly)",
                   color=C["bin"])

    y5 = 0.475
    bw2 = 0.55
    arrow(ax, 0.50, y4b - 0.008, 0.50, y5 + bh/2)
    box(ax, 0.50, y5, bw2, bh, C["bin"],
        "Coverage mapping  (bwa-mem2)",
        "Map Illumina reads → assembly  |  samtools sort  |  jgi_summarize_bam_contig_depths",
        fontsize=8)

    y6 = 0.410
    arrow(ax, 0.50, y5 - bh/2, 0.50, y6 + bh/2)
    box(ax, 0.50, y6, bw2, bh, C["bin"],
        "MetaBAT2  (per-source, multi-sample coverage)",
        "min_contig=2500  seed=1  threads=32\n(per-sample hybrid ×5 + co-assembly hybrid)",
        fontsize=8)
    param_note(ax, 0.92, y6 + bh/2 - 0.01,
               ["Bins from all","sources pooled","before QC"],
               fontsize=6)

    # ─── MAG QC section ────────────────────────────────────────────────────
    section_header(ax, 0.03, 0.383, "MAG QUALITY CONTROL & DEREPLICATION",
                   color=C["mag"])

    y7 = 0.350
    bw3 = 0.25
    boxes_qc = [
        (0.22, "CheckM2\n(completeness / contamination)",
         "threads=32  uniref100 KO diamond DB",
         "Filter: completeness ≥50%\ncontamination ≤10%"),
        (0.50, "dRep\n(genome dereplication)",
         "skani  ANI ≥95 %  -p 32",
         "276 dereplicated MAGs\n(species-level clusters)"),
        (0.78, "GTDB-Tk\n(taxonomy classification)",
         "classify_wf  threads=32\nrelease 220",
         "NCBI-compatible lineage\nassigned to each MAG"),
    ]
    prev_x = 0.22
    arrow(ax, 0.50, y6 - bh/2, 0.22, y7 + bh/2,
          connectionstyle="arc3,rad=0.15")
    arrow(ax, 0.50, y6 - bh/2, 0.50, y7 + bh/2)
    arrow(ax, 0.50, y6 - bh/2, 0.78, y7 + bh/2,
          connectionstyle="arc3,rad=-0.15")
    for xc, ttl, sub, note in boxes_qc:
        box(ax, xc, y7, bw3, bh * 1.2, C["mag"], ttl, sub, fontsize=7.5)
        param_note(ax, xc, y7 - bh * 0.6 - 0.02, [note], fontsize=6.5)

    arrow(ax, 0.22, y7 - bh * 0.6 - 0.01, 0.50, y7 - bh * 0.6 - 0.01)
    arrow(ax, 0.50, y7 - bh * 0.6 - 0.01, 0.50, 0.243)

    # ─── ANNOTATION section ────────────────────────────────────────────────
    section_header(ax, 0.03, 0.235, "MAG ANNOTATION  (276 dereplicated MAGs)",
                   color=C["ann"])

    y8 = 0.200
    bw4 = 0.26
    ann_boxes = [
        (0.22, "Bakta\n(gene prediction + function)",
         "Full DB (~70 GB)  threads=8×16 parallel\n--skip-plot\nOutputs: .gff3 .gbff .tsv .ffn .faa",
         C["ann"]),
        (0.50, "eggNOG-mapper\n(functional orthology)",
         "emapper.py  --itype genome\n--genepred prodigal  --sensmode fast\n--dmnd_iterate no  --dbmem\nthreads=8×16 parallel",
         C["ann"]),
        (0.78, "DRAM  (opt-in, disabled)\nMetabolism annotation",
         "~500 GB refdata  threads=32\nDRAM.py annotate\n[not run — opt-in only]",
         "#aaaaaa"),
    ]
    for xc, ttl, sub, col in ann_boxes:
        arrow(ax, 0.50, 0.243, xc, y8 + bh * 0.6,
              connectionstyle=f"arc3,rad={0.15 * (xc - 0.50):.2f}")
        box(ax, xc, y8, bw4, bh * 1.2, col, ttl, sub, fontsize=7, alpha=0.85)

    # ─── Outputs ──────────────────────────────────────────────────────────
    y9 = 0.090
    out_boxes = [
        (0.18, "ASV → MAG\n16S mapping", "minimap2 -ax sr\nconfidence tiers", C["abund"]),
        (0.40, "16S-bridge\nMAG abundance", "mag_abundance_by_day\n_union.csv  276×75 d", C["abund"]),
        (0.62, "Direct metagenomic\nMAG abundance", "bwa-mem2 + jgi_depth\nmag_abundance.csv 276×5", C["abund"]),
        (0.82, "Functional\nprofiles", "KEGG KO  COG\nCAZy  CARD", C["ann"]),
    ]
    for xc, ttl, sub, col in out_boxes:
        box(ax, xc, y9, 0.18, bh * 1.1, col, ttl, sub, fontsize=7.5)

    # Connect annotation to outputs
    for (xc, *_) in ann_boxes[:2]:
        arrow(ax, xc, y8 - bh * 0.6, 0.29, y9 + bh * 0.55,
              connectionstyle=f"arc3,rad={0.1 * np.sign(xc - 0.35):.2f}")
    arrow(ax, 0.50, y8 - bh * 0.6, 0.62, y9 + bh * 0.55,
          connectionstyle="arc3,rad=-0.15")
    arrow(ax, 0.50, y8 - bh * 0.6, 0.82, y9 + bh * 0.55,
          connectionstyle="arc3,rad=-0.25")

    # ─── Legend ───────────────────────────────────────────────────────────
    legend_items = [
        (C["qc"],    "QC / read prep"),
        (C["asm"],   "Assembly / polishing"),
        (C["bin"],   "Binning"),
        (C["mag"],   "MAG QC / taxonomy"),
        (C["ann"],   "Annotation"),
        (C["abund"], "Abundance"),
    ]
    lx, ly = 0.03, 0.055
    for col, lbl in legend_items:
        patch = FancyBboxPatch((lx, ly - 0.010), 0.018, 0.018,
                               boxstyle="round,pad=0.002",
                               facecolor=col, edgecolor="white", linewidth=1,
                               zorder=6)
        ax.add_patch(patch)
        ax.text(lx + 0.024, ly - 0.001, lbl, ha="left", va="center",
                fontsize=7, fontfamily=FONT, color=C["header"])
        lx += 0.155

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# PAGE 2 – 16S → MAG mapping  +  MAG abundance methods
# ─────────────────────────────────────────────────────────────────────────────

def make_page2(pdf):
    fig, axes = plt.subplots(1, 2, figsize=(17, 14),
                             gridspec_kw={"width_ratios": [1, 1]})
    page_background(fig)
    fig.suptitle("MAG Abundance Estimation — Two Approaches",
                 fontsize=14, fontfamily=FONT, fontweight="bold",
                 color=C["header"], y=0.98)

    # ── LEFT: 16S bridge ─────────────────────────────────────────────────
    ax = axes[0]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_facecolor(C["bg"])
    ax.text(0.5, 0.97, "Approach A — 16S Amplicon Bridge",
            ha="center", fontsize=12, fontfamily=FONT, fontweight="bold",
            color=C["abund"])
    ax.text(0.5, 0.945, "75 timepoints  (days 1–584)  |  indirect via DADA2 ASVs",
            ha="center", fontsize=8, fontfamily=FONT, color=C["note"])

    bw, bh = 0.68, 0.065
    steps_L = [
        (0.90, "Raw 16S amplicon reads\n(108 samples, paired-end)",
         "V4 region  515F/806R primers", C["qc"]),
        (0.800, "QIIME2 DADA2 denoising",
         "dada2 denoise-paired\n3950 ASVs × 75 samples\nrep_seqs_merged_dada2.qza", C["qc"]),
        (0.695, "QIIME2 relative abundance table",
         "table_merged_dada2.qza\nFeature table → rel. abundance per sample\ntable_rel_export.csv", C["qc"]),
        (0.585, "MAG 16S gene extraction",
         "Bakta .ffn files → grep 16S_rRNA\n343 genes from 239 / 276 MAGs\nmag_16s.fasta", C["mag"]),
        (0.475, "minimap2 alignment  (-ax sr)",
         "asvs.fasta  →  mag_16s.fasta\n3950 ASVs aligned (short-read preset)\nasv_vs_mag16s.sam", C["abund"]),
        (0.365, "Confidence tier assignment",
         "Species:  ≥97% id + ≥80% qcov  →  1339 ASVs\n"
         "Genus:    ≥95% id + ≥70% qcov  →   340 ASVs\n"
         "Family:   ≥90% id + ≥60% qcov  →   801 ASVs\n"
         "Weak:     below family          →  1318 ASVs\n"
         "No hit:                         →   152 ASVs", C["abund"]),
        (0.220, "ASV abundance → MAG transfer",
         "Each ASV's relative abundance assigned to\nbest-hit MAG at its confidence tier\nasv_to_mag_mapping.json", C["abund"]),
        (0.110, "mag_abundance_by_day_union.csv\nmag_abundance_by_day_intersection.csv",
         "276 MAGs × 75 days  |  values in %\nColumn sums ~99%  (1% = unmapped ASVs)\nUnion = all tiers  |  Intersection = high-conf only", C["abund"]),
    ]
    prev_y = None
    for y, ttl, sub, col in steps_L:
        if prev_y is not None:
            arrow(ax, 0.5, prev_y - bh / 2, 0.5, y + bh / 2)
        box(ax, 0.5, y, bw, bh, col, ttl, sub, fontsize=8)
        prev_y = y

    # Strengths / limitations
    ax.text(0.5, 0.025,
            "✓  75 timepoints   ✓  dense temporal resolution\n"
            "⚠  indirect (2-step error: DADA2 + alignment)   ⚠  37 MAGs have no 16S gene  →  always zero",
            ha="center", va="center", fontsize=7.5, fontfamily=FONT,
            color=C["header"],
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#eafaf1",
                      edgecolor=C["abund"], alpha=0.9))

    # ── RIGHT: direct metagenomics ────────────────────────────────────────
    ax = axes[1]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_facecolor(C["bg"])
    ax.text(0.5, 0.97, "Approach B — Direct Metagenomic Mapping",
            ha="center", fontsize=12, fontfamily=FONT, fontweight="bold",
            color="#16a085")
    ax.text(0.5, 0.945, "5 timepoints  (CAN 1–5)  |  direct from raw metagenome reads",
            ha="center", fontsize=8, fontfamily=FONT, color=C["note"])

    steps_R = [
        (0.88, "Illumina PE reads  (CAN 1–5)",
         "fastp-cleaned  R1 + R2 per sample\n~5 GB / sample", C["qc"]),
        (0.775, "Dereplicated MAG reference",
         "276 MAG FASTA files\n→ concatenated  mag/abundance/ref/mags.fasta\nbwa-mem2 index", C["mag"]),
        (0.665, "bwa-mem2 alignment  (per sample)",
         "bwa-mem2 mem -t 16\n→ samtools sort  →  per-sample .bam\n5 BAMs total", C["abund"]),
        (0.555, "jgi_summarize_bam_contig_depths",
         "All 5 BAMs → single depth.txt\nColumns: contigLen  totalAvgDepth\n  CAN_1  CAN_1-var  …  CAN_5  CAN_5-var", C["abund"]),
        (0.420, "compute_mag_abundance.py",
         "contig → MAG lookup from bin FASTA headers\n"
         "weighted mean depth =\n"
         "  Σ(depth_i × len_i) / Σ(len_i)  per MAG\n"
         "Fraction of binned bases:\n"
         "  pct = (depth × len_MAG) / Σ_all_MAGs(depth × len)\n"
         "RPKM = depth × 10⁹ / (150 bp × total_binned_bases)",
         "#16a085"),
        (0.265, "mag/abundance/mag_abundance.csv",
         "276 MAGs × 5 samples\nColumns per sample: mean_depth, pct, rpkm\n54.5% of contigs binned → denominator = all binned bases",
         "#16a085"),
    ]
    prev_y = None
    for y, ttl, sub, col in steps_R:
        if prev_y is not None:
            arrow(ax, 0.5, prev_y - bh / 2, 0.5, y + bh / 2)
        scale = 1.6 if "compute" in ttl else 1.0
        box(ax, 0.5, y, bw, bh * scale, col, ttl, sub, fontsize=7.5)
        prev_y = y - (bh * 0.3 if "compute" in ttl else 0)

    ax.text(0.5, 0.025,
            "✓  Direct from metagenome — no 16S intermediary\n"
            "✓  Captures all 276 MAGs including those without 16S genes\n"
            "⚠  Only 5 timepoints  |  ⚠  54.5% assembly binned (remainder unmeasured)\n"
            "→  Recommended use: validate 16S-bridge values at matched timepoints",
            ha="center", va="center", fontsize=7.5, fontfamily=FONT,
            color=C["header"],
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#e8f8f5",
                      edgecolor="#16a085", alpha=0.9))

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# PAGE 3 – Bakta vs alternative annotation tools
# ─────────────────────────────────────────────────────────────────────────────

def make_page3(pdf):
    fig, ax = plt.subplots(figsize=(17, 14))
    page_background(fig)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_facecolor(C["bg"])

    ax.text(0.5, 0.975, "Annotation Tool Selection — Bakta vs. Alternatives",
            ha="center", va="top", fontsize=15, fontfamily=FONT,
            fontweight="bold", color=C["header"])
    ax.text(0.5, 0.950, "Justification for using Bakta (full DB) as the primary annotation tool for 276 dereplicated MAGs",
            ha="center", va="top", fontsize=9, fontfamily=FONT, color=C["note"])

    # ── Comparison table ──────────────────────────────────────────────────
    tools = ["Bakta\n(chosen)", "Prokka", "RAST\n(SEED)", "DRAM", "Prodigal\n(ORF only)"]
    criteria = [
        "Database",
        "Reference currency",
        "16S rRNA\nrecovery",
        "Reproducibility\n/ local run",
        "Runtime\n(×276 MAGs)",
        "Output formats",
        "Functional\ndepth",
        "MAG / draft\ngenome support",
    ]
    cell_data = {
        # (tool_idx, criterion_idx): (text, color)  color=None → default
        # Bakta
        (0, 0): ("RefSeq NR, UniProt SwissProt,\nISfinder, CARD, REBASE, AMRFinder,\nSignalP", None),
        (0, 1): ("Updated to RefSeq 2024\n(curated, versioned)", "good"),
        (0, 2): ("✓ Explicit 16S/23S/5S\nrRNA annotation\n→ used for ASV–MAG mapping", "good"),
        (0, 3): ("✓ Fully local,\ndocker/conda,\nversion-pinned", "good"),
        (0, 4): ("~2 h total\n(16 parallel × 8 t)", "good"),
        (0, 5): (".gff3 .gbff .tsv\n.ffn .faa .log\n(INSDC-compliant)", "good"),
        (0, 6): ("Gene function, EC,\nCOG, GO, eggNOG,\nARG, IS elements", "good"),
        (0, 7): ("✓ Designed for MAGs\n& draft genomes", "good"),

        # Prokka
        (1, 0): ("UniProt + HAMAP\n(smaller, older)", None),
        (1, 1): ("Databases last updated\n~2021; lags RefSeq", "warn"),
        (1, 2): ("✓ rRNA annotated\nbut less curated", "ok"),
        (1, 3): ("✓ Local / conda", "good"),
        (1, 4): ("Similar to Bakta", "ok"),
        (1, 5): (".gff .gbk .faa .ffn\n.tsv .txt", "ok"),
        (1, 6): ("Gene function + COG\n(no CARD/IS/AMR\nby default)", "warn"),
        (1, 7): ("✓ Widely used for MAGs\nbut showing age", "ok"),

        # RAST
        (2, 0): ("SEED subsystems\n(curated, large)", None),
        (2, 1): ("SEED updated, but\nclosed-access", "warn"),
        (2, 2): ("✓ Included in\nRASTtk pipeline", "ok"),
        (2, 3): ("✗ Web-based API\nor RASTtk server;\nqueue delays,\nno versioning", "bad"),
        (2, 4): ("Web queue: days\nRASTtk local: OK\nbut setup complex", "warn"),
        (2, 5): (".gff .gbk .tsv\nSEED-specific IDs", "ok"),
        (2, 6): ("Subsystem-based;\ngood metabolic\ncontext", "ok"),
        (2, 7): ("✓ Works; not\noptimised for MAGs", "ok"),

        # DRAM
        (3, 0): ("KEGG, UniRef90,\nPfam, CARD, dbCAN,\nVOG (~500 GB)", None),
        (3, 1): ("Comprehensive but\nrequires large\ndisk & RAM", "ok"),
        (3, 2): ("✓ Annotates 16S\nif provided", "ok"),
        (3, 3): ("✓ Local conda\nbut heavy setup\n(500 GB refdata)", "warn"),
        (3, 4): ("Very slow:\nhours/MAG\n276 MAGs = weeks", "bad"),
        (3, 5): (".tsv distill report\n+ genome summaries", "ok"),
        (3, 6): ("Best metabolism\ncontext; KEGG\npathway maps", "good"),
        (3, 7): ("✓ Designed for MAGs;\ncomplex distillation\nstep required", "ok"),

        # Prodigal
        (4, 0): ("None — ORF\nprediction only", None),
        (4, 1): ("N/A", None),
        (4, 2): ("✗ No rRNA;\nexternal tool\nneeded", "bad"),
        (4, 3): ("✓ Ultrafast,\nno DB needed", "good"),
        (4, 4): ("< 1 min / MAG", "good"),
        (4, 5): (".gff .faa .fna", "ok"),
        (4, 6): ("ORF boundaries\nonly; no function", "bad"),
        (4, 7): ("✓ Core of most\ntools (Bakta\nuses it internally)", "ok"),
    }

    color_map = {
        "good": "#d5f5e3",
        "ok":   "#fef9e7",
        "warn": "#fde8d8",
        "bad":  "#fadbd8",
        None:   "#ffffff",
    }

    # Table geometry
    col_w = [0.16, 0.12, 0.12, 0.12, 0.12, 0.12]   # crit-col + 5 tool cols
    col_x = [0.01]
    for w in col_w[1:]:
        col_x.append(col_x[-1] + col_w[0] + (w if col_x[-1] == 0.01 else w))
    # recompute properly
    col_x = [0.01]
    for i in range(1, len(tools) + 1):
        prev = col_x[-1]
        w = col_w[0] if i == 0 else (0.84 / len(tools))
        col_x.append(prev + (col_w[0] if i == 1 else 0.84 / len(tools)))
    # simple grid
    crit_col_x = 0.01
    tool_xs = [0.01 + col_w[0] + (0.84 / len(tools)) * (i + 0.5) for i in range(len(tools))]
    row_h = 0.080
    table_top = 0.905
    table_left = 0.01
    table_width = 0.98
    crit_col_w = col_w[0]

    # Header row — tool names
    y_hdr = table_top
    rect = FancyBboxPatch((table_left, y_hdr - 0.030), table_width, 0.030,
                          boxstyle="round,pad=0.002",
                          facecolor=C["header"], edgecolor="none", zorder=3)
    ax.add_patch(rect)
    for xi, tool in zip(tool_xs, tools):
        weight = "bold" if "Bakta" in tool else "normal"
        ax.text(xi, y_hdr - 0.015, tool,
                ha="center", va="center", fontsize=9, fontfamily=FONT,
                fontweight=weight, color="white", zorder=4)

    # Rows
    for ri, crit in enumerate(criteria):
        y_top = table_top - 0.030 - ri * row_h
        y_cen = y_top - row_h / 2
        bg = "#f2f3f4" if ri % 2 == 0 else "#ffffff"
        rect = FancyBboxPatch((table_left, y_top - row_h), table_width, row_h,
                              boxstyle="round,pad=0.001",
                              facecolor=bg, edgecolor="#ddd", linewidth=0.5, zorder=2)
        ax.add_patch(rect)
        # Criterion label
        ax.text(crit_col_x + crit_col_w / 2, y_cen, crit,
                ha="center", va="center", fontsize=8, fontfamily=FONT,
                fontweight="bold", color=C["header"], zorder=4,
                multialignment="center")
        # Cell data
        for ci, xi in enumerate(tool_xs):
            key = (ci, ri)
            if key in cell_data:
                txt, rating = cell_data[key]
                cell_col = color_map[rating]
                cell_w = 0.84 / len(tools) - 0.005
                cell_x = xi - cell_w / 2
                if cell_col != "#ffffff":
                    crect = FancyBboxPatch(
                        (cell_x, y_top - row_h + 0.004), cell_w, row_h - 0.008,
                        boxstyle="round,pad=0.003",
                        facecolor=cell_col, edgecolor="none", zorder=3, alpha=0.85,
                    )
                    ax.add_patch(crect)
                ax.text(xi, y_cen, txt,
                        ha="center", va="center", fontsize=6.8,
                        fontfamily=FONT, color=C["header"],
                        zorder=5, multialignment="center")

    # ── Why Bakta text ─────────────────────────────────────────────────────
    y_text = table_top - 0.030 - len(criteria) * row_h - 0.025
    rationale = (
        "Why Bakta was chosen for this pipeline\n"
        "─────────────────────────────────────────────────────────────────────────────────────────────────\n"
        "1. Explicit rRNA annotation (16S / 23S / 5S):  Bakta outputs .ffn files containing all non-coding RNA sequences, including annotated\n"
        "   16S rRNA genes.  These were used directly to build the MAG–16S reference (mag_16s.fasta) for the ASV → MAG abundance bridge.\n"
        "   Prokka also annotates rRNA, but Bakta's 16S loci carry RefSeq-backed gene IDs that make downstream matching more reliable.\n\n"
        "2. Currency of reference databases:  Bakta ships versioned snapshots of RefSeq NR, UniProt Swiss-Prot, ISfinder, CARD, REBASE,\n"
        "   and AMRFinderPlus.  Prokka's bundled databases were last comprehensively updated in 2021.  For a 2026 wastewater microbiome\n"
        "   project (dominated by EBPR taxa with many 'hypothetical' proteins), current databases significantly improve hit rates.\n\n"
        "3. Local, reproducible, version-pinned:  RAST requires either a web API (introducing queue latency and non-reproducibility)\n"
        "   or a self-hosted RASTtk instance.  Bakta runs offline via conda with an explicit DB version — essential for a shared HPC\n"
        "   environment where connectivity and API stability cannot be guaranteed.\n\n"
        "4. DRAM not run (opted out):  DRAM offers the richest metabolic reconstruction but requires ~500 GB of reference data, takes\n"
        "   hours per MAG (weeks for 276 MAGs serial), and is primarily useful for metabolic pathway distillation — a downstream\n"
        "   analysis step.  The Bakta + eggNOG combination (KEGG KO, COG, CAZy, GO, EC numbers) provides equivalent functional\n"
        "   coverage for the community ecology questions (GAO/PAO enrichment, nutrient cycling) at a fraction of the compute cost.\n\n"
        "5. INSDC-compliant output:  Bakta produces .gbff files that can be deposited directly to NCBI/ENA.  This is required for\n"
        "   the KBase MAG + genome import workflow already configured in kbase_mag_import.csv and kbase_genome_import.csv."
    )
    ax.text(0.01, y_text, rationale,
            ha="left", va="top", fontsize=8, fontfamily="DejaVu Sans Mono",
            color=C["header"], zorder=5, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#fdfefe",
                      edgecolor=C["ann"], alpha=0.92))

    # Colour legend for table
    lx = 0.01
    ly = 0.025
    for col, lbl in [(color_map["good"], "Better than alternatives"),
                     (color_map["ok"],   "Adequate"),
                     (color_map["warn"], "Some limitation"),
                     (color_map["bad"],  "Significant limitation")]:
        patch = FancyBboxPatch((lx, ly - 0.010), 0.022, 0.016,
                               boxstyle="round,pad=0.002",
                               facecolor=col, edgecolor="#aaa", linewidth=0.8, zorder=6)
        ax.add_patch(patch)
        ax.text(lx + 0.027, ly - 0.002, lbl, ha="left", va="center",
                fontsize=7.5, fontfamily=FONT, color=C["header"])
        lx += 0.230

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    plt.rcParams.update({
        "pdf.fonttype": 42,    # embed TrueType fonts
        "ps.fonttype":  42,
        "figure.dpi":   150,
    })
    with PdfPages(OUT) as pdf:
        # PDF metadata
        d = pdf.infodict()
        d["Title"]   = "EmilyKin Metagenomics Pipeline Diagram"
        d["Author"]  = "Andrew Freiburger"
        d["Subject"] = "Hybrid ONT+Illumina metagenomics: assembly → MAGs → annotation → abundance"
        d["Keywords"] = "metagenomics Flye Medaka MetaBAT2 dRep GTDB-Tk Bakta eggNOG MAG abundance"

        print("Rendering page 1 (pipeline)...")
        make_page1(pdf)
        print("Rendering page 2 (abundance methods)...")
        make_page2(pdf)
        print("Rendering page 3 (Bakta justification)...")
        make_page3(pdf)

    print(f"\nDone → {OUT}")


if __name__ == "__main__":
    main()
