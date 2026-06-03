#!/usr/bin/env python3
"""
Map raw 16S amplicon reads directly to MAG-extracted 16S genes (Approach 2).

For each of the 75 study samples, maps R1+R2 reads to mag_16s.fasta using
minimap2 (-ax sr, primary alignments only), counts reads per gene with
samtools idxstats, aggregates by MAG, and normalises to relative abundance.

Output: mag_16s_abundance.csv  (276 MAGs × 75 samples, values in %)

Run from the EmilyKin directory:
    python map_16s_to_mags.py
"""
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
RAWREADS   = Path("/scratch1/afreiburger/emilykin/raw/16s/rawreads")
MAG16S     = Path("/scratch1/afreiburger/emilykin/asv_mag_mapping/mag_16s.fasta")
SAMPLE_DAYS = Path("sample_days.json")
OUT_CSV    = Path("mag_16s_abundance.csv")
WORK_DIR   = Path("mag_16s_mapping_tmp")

# assembly conda env minimap2 + samtools
MINIMAP2   = "/scratch/afreiburger/metag-hybrid/env/miniforge/envs/assembly/bin/minimap2"
SAMTOOLS   = "/scratch/afreiburger/metag-hybrid/env/miniforge/envs/assembly/bin/samtools"

THREADS_PER_JOB = 4   # minimap2 -t per sample
MAX_PARALLEL    = 16  # concurrent sample jobs


# ── helpers ──────────────────────────────────────────────────────────────────

def parse_gene_to_mag(fasta: Path) -> dict[str, str]:
    """Return {gene_id: mag_name} by parsing headers like '>MAG__LOCUS_TAG ...'"""
    gene_to_mag = {}
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            # header format: >CAN_1_bin.148__GDJGFH_00714 [optional description]
            header = line[1:].split()[0]
            if "__" in header:
                mag, gene = header.split("__", 1)
            else:
                mag, gene = header, header
            gene_to_mag[header] = mag
    return gene_to_mag


def build_sample_map(rawreads: Path, study_samples: set) -> dict[str, tuple]:
    """Return {sample_id: (r1_path, r2_path)} for all study samples."""
    r1, r2 = {}, {}
    for f in rawreads.iterdir():
        m = re.search(r'EK[-_](B\w+)_S', f.name)
        if not m:
            continue
        sid = m.group(1)
        if sid not in study_samples:
            continue
        if "_R1_" in f.name:
            r1[sid] = f
        elif "_R2_" in f.name:
            r2[sid] = f
    paired = {s: (r1[s], r2[s]) for s in study_samples if s in r1 and s in r2}
    return paired


def map_sample(sample_id: str, r1: Path, r2: Path,
               ref: Path, work: Path) -> dict[str, int]:
    """
    Map one sample with minimap2 and return {gene_id: read_count}.
    Counts primary, mapped, first-in-pair reads only (avoids double-counting PE).
    """
    bam = work / f"{sample_id}.bam"
    log = work / f"{sample_id}.log"

    with open(log, "w") as lf:
        # minimap2: short-read preset, no secondary alignments, pipe to samtools
        # -F 2308 = exclude: unmapped(4) + not-primary(256) + supplementary(2048)
        # -f 1    = include: paired (retain only first-in-pair for counting)
        # Using -f 64 (first-in-pair) to count each fragment once
        cmd = (
            f"{MINIMAP2} -ax sr --secondary=no -t {THREADS_PER_JOB} "
            f"{ref} {r1} {r2} 2>>{log} "
            f"| {SAMTOOLS} view -bF 2308 -f 64 "
            f"| {SAMTOOLS} sort -@ 2 -o {bam} -"
        )
        ret = subprocess.run(cmd, shell=True, stderr=lf)
        if ret.returncode != 0:
            raise RuntimeError(f"{sample_id}: minimap2/samtools failed (rc={ret.returncode})")

        subprocess.run([SAMTOOLS, "index", str(bam)], check=True, stderr=lf)

        # idxstats: ref_name, ref_len, mapped, unmapped — one line per reference
        idxstats = subprocess.run(
            [SAMTOOLS, "idxstats", str(bam)],
            capture_output=True, text=True, check=True
        )

    counts = {}
    for line in idxstats.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        gene_id, _ref_len, mapped, *_ = parts
        if gene_id == "*":
            continue
        counts[gene_id] = int(mapped)

    # clean up BAM to save space (keep log)
    bam.unlink(missing_ok=True)
    Path(str(bam) + ".bai").unlink(missing_ok=True)

    return counts


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    WORK_DIR.mkdir(exist_ok=True)

    print("Parsing MAG 16S gene map...", flush=True)
    gene_to_mag = parse_gene_to_mag(MAG16S)
    all_mags = sorted(set(gene_to_mag.values()))
    print(f"  {len(gene_to_mag)} genes → {len(all_mags)} MAGs", flush=True)

    print("Loading sample metadata...", flush=True)
    sd = json.load(open(SAMPLE_DAYS))
    study_samples = set(sd.keys())
    sample_map = build_sample_map(RAWREADS, study_samples)
    print(f"  {len(sample_map)} / {len(study_samples)} samples have R1+R2", flush=True)
    if len(sample_map) < len(study_samples):
        missing = study_samples - set(sample_map)
        print(f"  WARNING: missing FASTQs for: {missing}", file=sys.stderr)

    # results[sample][mag] = relative_abundance (%)
    results: dict[str, dict[str, float]] = {}

    print(f"\nMapping {len(sample_map)} samples (max {MAX_PARALLEL} parallel, "
          f"{THREADS_PER_JOB} threads each)...", flush=True)

    def run_one(args):
        sid, r1, r2 = args
        counts = map_sample(sid, r1, r2, MAG16S, WORK_DIR)
        # aggregate counts by MAG
        mag_counts: dict[str, int] = defaultdict(int)
        for gene_id, cnt in counts.items():
            mag = gene_to_mag.get(gene_id, gene_id)
            mag_counts[mag] += cnt
        total = sum(mag_counts.values())
        if total == 0:
            return sid, {mag: 0.0 for mag in all_mags}
        return sid, {mag: mag_counts.get(mag, 0) / total * 100 for mag in all_mags}

    jobs = [(sid, r1, r2) for sid, (r1, r2) in sample_map.items()]
    done = 0
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as pool:
        futures = {pool.submit(run_one, j): j[0] for j in jobs}
        for fut in as_completed(futures):
            sid = futures[fut]
            try:
                sid_out, abund = fut.result()
                results[sid_out] = abund
                done += 1
                total_mapped = sum(abund.values())
                print(f"  [{done}/{len(jobs)}] {sid_out}: "
                      f"{total_mapped:.1f}% reads mapped to MAG 16S genes", flush=True)
            except Exception as exc:
                print(f"  ERROR {sid}: {exc}", file=sys.stderr)

    # Write output CSV: rows=MAGs, cols=samples sorted by day
    print(f"\nWriting {OUT_CSV} ...", flush=True)
    samples_sorted = sorted(results.keys(), key=lambda s: sd[s])

    with open(OUT_CSV, "w") as f:
        day_header = "\t".join(f"day_{sd[s]}" for s in samples_sorted)
        f.write(f"MAG\t{day_header}\n")
        for mag in all_mags:
            vals = "\t".join(f"{results[s].get(mag, 0.0):.6f}" for s in samples_sorted)
            f.write(f"{mag}\t{vals}\n")

    print(f"Done → {OUT_CSV}  ({len(all_mags)} MAGs × {len(results)} samples)")
    print(f"Logs in {WORK_DIR}/")


if __name__ == "__main__":
    main()
