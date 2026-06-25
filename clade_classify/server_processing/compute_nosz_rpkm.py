#!/usr/bin/env python3
"""
SERVER-SIDE nosZ RPKM + depth-based copy number, split by Clade I/II/III, for ALL
nosZ-bearing bins.  Run on poplar (where the sample BAMs + samtools live).

Why this runs on the server, not the workstation: it needs the per-sample read
alignments (aligned_<sample>_sorted.bam, mapped to the dereplicated-MAG combined
reference) and samtools — the same inputs the original gene-abundance step used
(gene_ab_figure/05_gene_copy_number.py, 03_count_and_rpkm.py).

Dependencies: python3 (standard library only) + samtools. No pip installs.

Inputs (in this directory unless you edit CONFIG):
  allbins_nosz_loci_clade.bed   <- produced locally; nosZ loci with BAM-style contig
                                   names, coordinates, MAG, and Clade I/II/III call.
  aligned_<sample>_sorted.bam (+ .bai) for each sample.

Outputs (OUTDIR):
  nosz_locus_rpkm.tsv              per-locus x per-sample: reads, RPKM, copy number
  nosz_rpkm_per_sample_claded.tsv  per-bin x per-sample, long format
                                   (ko_group = nosZ_cladeI/II/III) — same layout as
                                   gene_ab_figure/data/gene_rpkm_per_sample_claded.tsv
  nosz_copynumber_per_sample.tsv   per-bin x per-sample depth-based copy number by clade
  sample_totals.tsv                total primary-mapped reads per sample

Run:  python3 compute_nosz_rpkm.py
"""
import csv, os, subprocess
from collections import defaultdict

# ============================ CONFIG — edit to match the server ============================
SAMTOOLS = "/scratch1/afreiburger/emilykin/processed/.snakemake_envs/88fdb48d4d745c55ec2cd90b407de422_/bin/samtools"
BAM_DIR  = "/scratch1/afreiburger/emilykin/gene_ab_figure/data"   # holds aligned_<sample>_sorted.bam
BAM_TMPL = "aligned_{sample}_sorted.bam"
SAMPLES  = ["CAN_1", "CAN_2", "CAN_3", "CAN_4", "CAN_5"]
BED      = "allbins_nosz_loci_clade.bed"
OUTDIR   = "nosz_rpkm_out"
FLAGS    = "2308"   # samtools -F 2308 = drop unmapped(4)+secondary(256)+supplementary(2048): primary mapped only
# ==========================================================================================

def sh(cmd):
    return subprocess.run([str(x) for x in cmd], capture_output=True, text=True, check=True).stdout

def main():
    if not os.path.exists(SAMTOOLS):
        raise SystemExit(f"samtools not found at {SAMTOOLS} — edit CONFIG.SAMTOOLS")
    os.makedirs(OUTDIR, exist_ok=True)

    loci = []
    with open(BED) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            loci.append({"contig": p[0], "start0": int(p[1]), "end": int(p[2]),
                         "locus_id": p[3], "mag": p[4], "clade": p[5]})
    print(f"{len(loci)} nosZ loci from {BED}")

    # ---- per-sample idxstats: contig lengths, per-bin mapped reads, total mapped ----
    contig_len = {}
    genome_reads = defaultdict(lambda: defaultdict(int))   # bin_prefix -> sample -> reads
    total = defaultdict(int)                                # sample -> total mapped
    bam_contigs = set()
    for s in SAMPLES:
        bam = os.path.join(BAM_DIR, BAM_TMPL.format(sample=s))
        if not os.path.exists(bam):
            raise SystemExit(f"missing BAM: {bam} — edit CONFIG.BAM_DIR/BAM_TMPL")
        if not (os.path.exists(bam + ".bai") or os.path.exists(bam[:-4] + ".bai")):
            print(f"  indexing {bam} ..."); sh([SAMTOOLS, "index", bam])
        for line in sh([SAMTOOLS, "idxstats", bam]).splitlines():
            c, ln, mapped, _ = line.split("\t")
            if c == "*":
                continue
            bam_contigs.add(c); contig_len[c] = int(ln)
            genome_reads[c.split("::")[0]][s] += int(mapped)
            total[s] += int(mapped)
        print(f"  {s}: total primary-mapped reads = {total[s]:,}")
    genome_len = defaultdict(int)
    for c, ln in contig_len.items():
        genome_len[c.split("::")[0]] += ln

    # ---- naming sanity check ----
    missing = sorted({L["contig"] for L in loci if L["contig"] not in bam_contigs})
    if missing:
        print(f"\n*** WARNING: {len(missing)}/{len(loci)} loci contigs are NOT in the BAM "
              f"(they will count 0). Example: {missing[:3]}")
        print("    The BED uses combined_ref.fa naming  <mag_with_underscores>::<contig>.")
        print("    If your BAM headers differ, fix CONFIG or the BED contig column.\n")
    else:
        print("naming check OK: all loci contigs found in the BAMs.\n")

    # ---- per-locus read counts -> RPKM + copy number ----
    rows = []
    for L in loci:
        region = f"{L['contig']}:{L['start0'] + 1}-{L['end']}"
        glen = L["end"] - L["start0"]; binp = L["contig"].split("::")[0]
        rec = {**L, "gene_len": glen}
        for s in SAMPLES:
            bam = os.path.join(BAM_DIR, BAM_TMPL.format(sample=s))
            n = 0
            if L["contig"] in bam_contigs:
                n = int(sh([SAMTOOLS, "view", "-c", "-F", FLAGS, bam, region]).strip())
            rpkm = n / (glen / 1000.0) / (total[s] / 1e6) if (total[s] and glen) else 0.0
            gr, gln = genome_reads[binp][s], genome_len[binp]
            cn = (n / glen) / (gr / gln) if (gr and gln and glen) else 0.0
            rec[f"{s}_reads"], rec[f"{s}_rpkm"], rec[f"{s}_copynum"] = n, f"{rpkm:.4f}", f"{cn:.4f}"
        rows.append(rec)

    cols = (["locus_id", "mag", "clade", "contig", "start0", "end", "gene_len"]
            + [f"{s}_{x}" for s in SAMPLES for x in ("reads", "rpkm", "copynum")])
    with open(os.path.join(OUTDIR, "nosz_locus_rpkm.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    # ---- per-bin, split by clade (sum a bin's same-clade loci) ----
    bins = sorted({r["mag"] for r in rows})
    rpkm_bc = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))  # clade->mag->sample->rpkm
    cn_bc   = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    for r in rows:
        for s in SAMPLES:
            rpkm_bc[r["clade"]][r["mag"]][s] += float(r[f"{s}_rpkm"])
            cn_bc[r["clade"]][r["mag"]][s]   += float(r[f"{s}_copynum"])
    def write_claded(path, data):
        with open(path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t"); w.writerow(["ko_group", "mag"] + SAMPLES)
            for clade in ("I", "II", "III"):
                for mag in bins:
                    if mag in data[clade]:
                        w.writerow([f"nosZ_clade{clade}", mag]
                                   + [f"{data[clade][mag][s]:.4f}" for s in SAMPLES])
    write_claded(os.path.join(OUTDIR, "nosz_rpkm_per_sample_claded.tsv"), rpkm_bc)
    write_claded(os.path.join(OUTDIR, "nosz_copynumber_per_sample.tsv"), cn_bc)

    with open(os.path.join(OUTDIR, "sample_totals.tsv"), "w") as f:
        f.write("sample\ttotal_primary_mapped_reads\n")
        for s in SAMPLES:
            f.write(f"{s}\t{total[s]}\n")

    print(f"\nDONE -> {OUTDIR}/")
    print("  nosz_locus_rpkm.tsv  nosz_rpkm_per_sample_claded.tsv  "
          "nosz_copynumber_per_sample.tsv  sample_totals.tsv")

if __name__ == "__main__":
    main()
