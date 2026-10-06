#!/usr/bin/env python3
"""Stream samtools depth -aa stdin; compute primary mapping metrics from BAM."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("bam", "output", "platform", "seed", "fraction"):
        parser.add_argument("--" + field, required=True)
    args = parser.parse_args()
    positions = bases_at_depth = covered1 = covered5 = covered10 = 0
    for line in sys.stdin:
        columns = line.rstrip().split("\t")
        if len(columns) < 3:
            raise ValueError("Expected samtools depth chromosome/position/depth columns")
        depth = int(columns[2])
        positions += 1
        bases_at_depth += depth
        covered1 += depth >= 1
        covered5 += depth >= 5
        covered10 += depth >= 10
    if not positions:
        raise ValueError("samtools depth produced no positions; verify reference header and -aa flag")
    reads = mapped = nominal_bases = mapq_sum = 0
    error_fraction_sum = 0.0
    error_fraction_count = 0
    process = subprocess.Popen(["samtools", "view", "-F", "0x900", args.bam], stdout=subprocess.PIPE, text=True)
    for line in process.stdout:
        fields = line.rstrip().split("\t")
        reads += 1
        if fields[9] != "*":
            nominal_bases += len(fields[9])
        if not int(fields[1]) & 4:
            mapped += 1
            mapq_sum += int(fields[4])
            aligned_query = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]) if op in "MI=X")
            nm = next((int(t[5:]) for t in fields[11:] if t.startswith("NM:i:")), None)
            if nm is not None:
                error_fraction_sum += nm / max(1, aligned_query)
                error_fraction_count += 1
    if process.wait() != 0:
        raise RuntimeError("samtools view failed")
    report = {"platform": args.platform, "seed": int(args.seed), "fraction": float(args.fraction),
              "reference_positions": positions, "reads": reads, "mapped_reads": mapped,
              "nominal_depth_x": nominal_bases / positions, "mapped_mean_depth_x": bases_at_depth / positions,
              "mapping_rate": mapped / max(reads, 1), "mean_mapq_mapped": mapq_sum / max(mapped, 1),
              "mean_alignment_edit_distance_per_aligned_query_base": error_fraction_sum / error_fraction_count if error_fraction_count else None,
              "breadth_1x": covered1 / positions, "breadth_5x": covered5 / positions, "breadth_10x": covered10 / positions,
              "samtools_version": subprocess.check_output(["samtools", "--version"], text=True).splitlines()[0],
              "coverage_policy": "Primary alignments; SAM CIGAR matching positions, gaps excluded; -Q20 MAPQ by default; paired overlaps counted; duplicates not marked. samtools depth excludes preflagged QCFAIL/DUP; mapping-rate counts all primary reads. NM/(M+I+=+X) is an edit-distance/query-base proxy, not sequencing-error truth."}
    Path(args.output).write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
