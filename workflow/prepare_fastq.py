#!/usr/bin/env python3
"""ONT-specific streaming Phred+33 QC; no short-read adapter/tail assumptions."""
import argparse
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import QC, fastq_records, filter_read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    before, after, rejected = QC(), QC(), {}
    with gzip.open(args.output, "wt", encoding="ascii") as target:
        for read in fastq_records(args.input):
            before.add(read)
            filtered, reason = filter_read(read, "nanopore")
            if filtered is None:
                rejected[reason] = rejected.get(reason, 0) + 1
            else:
                after.add(filtered)
                target.write(f"@{filtered.name}\n{filtered.sequence}\n+\n{filtered.quality}\n")
    audit = {"platform": "nanopore", "thresholds": {"minimum_length": 500, "minimum_mean_phred": 7},
             "raw_qc": before.report(), "filtered_qc": after.report(), "rejections": rejected,
             "caveat": "Historical ONT quality scores can be poorly calibrated; this permissive threshold is not proof of true base accuracy."}
    Path(args.audit).write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
