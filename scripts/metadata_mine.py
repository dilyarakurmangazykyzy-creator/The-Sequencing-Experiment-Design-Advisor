#!/usr/bin/env python3
"""Stream ENA/SRA TSV metadata, audit missingness and aggregate archived designs.

This counts runs, not studies, biological samples or instrument sales. Archive
submission year is NOT the date an experiment was performed. No deduplication
between ENA/SRA mirrors is attempted: merge by run accession before combining.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import process_peak_memory_mb

ALIASES = {
    "run_accession": ("run_accession", "Run", "run"),
    "platform": ("instrument_platform", "Platform", "platform"),
    "model": ("instrument_model", "Model", "model"),
    "strategy": ("library_strategy", "LibraryStrategy", "strategy"),
    "layout": ("library_layout", "LibraryLayout", "layout"),
    "date": ("first_public", "published", "LoadDate", "ReleaseDate"),
    "bases": ("base_count", "bases", "Bases"),
    "reads": ("read_count", "reads", "Reads"),
    "spots": ("spot_count", "spots", "Spots"),
    "organism": ("scientific_name", "ScientificName", "organism"),
    "bioproject": ("study_accession", "BioProject", "Study"),
    "biosample": ("sample_accession", "BioSample", "Sample"),
}

FETCH_FIELDS = ["run_accession", "study_accession", "sample_accession", "scientific_name",
                "instrument_platform", "instrument_model", "library_strategy", "library_layout",
                "first_public", "read_count", "base_count"]


def fetch_metadata(destination, query, limit=1000000):
    """Stream the official ENA Portal API response to disk with provenance.

    ENA mirrors INSDC/SRA raw runs; this is not a second independent database.
    A bounded archive response is a query cohort, not a probability sample.
    """
    if limit < 0:
        raise ValueError("limit must be nonnegative; ENA limit=0 means all")
    parameters = {"result": "read_run", "query": query, "fields": ",".join(FETCH_FIELDS),
                  "format": "tsv", "limit": limit}
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urlencode(parameters)
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    digest = hashlib.sha256()
    received = 0
    with urlopen(Request(url, headers={"User-Agent": "AITU-Sequencing-Advisor/1.0"}), timeout=180) as response, partial.open("wb") as output:
        while block := response.read(1024 * 1024):
            output.write(block)
            digest.update(block)
            received += len(block)
    with partial.open(encoding="utf-8-sig") as check:
        header = check.readline()
    if "run_accession" not in header.split("\t"):
        raise ValueError("ENA response was not the requested metadata TSV; inspect .part response")
    partial.replace(target)
    provenance = {"url": url, "parameters": parameters, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                  "bytes": received, "sha256": digest.hexdigest(),
                  "scope": "Bounded ENA Portal API query response; actual rows are measured during aggregation.",
                  "documentation": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html"}
    target.with_suffix(target.suffix + ".provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return provenance


def get(row, field):
    for alias in ALIASES[field]:
        value = row.get(alias, "").strip()
        if value and value.lower() not in ("na", "n/a", "null", "none"):
            return value
    return ""


def positive_number(value):
    try:
        number = int(value.replace(",", ""))
        return number if number >= 0 else None
    except (ValueError, AttributeError):
        return None


def aggregate(input_path, output_dir, max_records=None, quiet=False):
    start = time.perf_counter()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    missing, invalid = Counter(), Counter()
    groups = defaultdict(lambda: {"runs": 0, "bases": 0, "reads": 0, "spots": 0,
                                  "runs_with_bases": 0, "runs_with_reads": 0, "runs_with_spots": 0,
                                  "bases_per_archive_unit_sum": 0.0, "runs_with_count_unit_ratio": 0,
                                  "complete_record_bases": 0, "complete_record_reads": 0})
    platforms, strategies, layouts, organisms = Counter(), Counter(), Counter(), Counter()
    examples = []
    rows = 0
    path = Path(input_path)
    opener = gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix == ".gz" else path.open(encoding="utf-8-sig")
    with opener as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Metadata TSV has no header")
        columns = reader.fieldnames
        schema_unavailable = [field for field, aliases in ALIASES.items() if not any(alias in columns for alias in aliases)]
        for row in reader:
            if max_records is not None and rows >= max_records:
                break
            rows += 1
            values = {field: get(row, field) for field in ALIASES}
            for field, value in values.items():
                missing[field] += not bool(value)
            if not re.fullmatch(r"[SED]RR\d+", values["run_accession"]):
                invalid["run_accession_format"] += 1
            year = values["date"][:4] if re.match(r"(?:19|20)\d{2}", values["date"]) else "unknown"
            if values["date"] and year == "unknown":
                invalid["submission_date"] += 1
            platform = values["platform"].upper() or "unknown"
            strategy = values["strategy"].upper() or "unknown"
            layout = values["layout"].upper() or "unknown"
            group = groups[(year, platform, strategy, layout)]
            group["runs"] += 1
            for field in ("bases", "reads", "spots"):
                n = positive_number(values[field])
                if n is not None:
                    group[field] += n
                    group["runs_with_" + field] += 1
                elif values[field]:
                    invalid[field + "_nonnegative_integer"] += 1
            bases, reads = positive_number(values["bases"]), positive_number(values["reads"])
            if bases is not None and reads:
                group["bases_per_archive_unit_sum"] += bases / reads
                group["runs_with_count_unit_ratio"] += 1
                group["complete_record_bases"] += bases
                group["complete_record_reads"] += reads
            platforms[platform] += 1
            strategies[strategy] += 1
            layouts[layout] += 1
            if values["organism"]:
                organisms[values["organism"]] += 1
            if len(examples) < 25:
                examples.append({"run_accession": values["run_accession"], "bioproject": values["bioproject"],
                                 "biosample": values["biosample"], "organism": values["organism"],
                                 "platform": platform, "strategy": strategy})
    design_rows = []
    for (year, platform, strategy, layout), counts in sorted(groups.items()):
        length_sum = counts.pop("bases_per_archive_unit_sum")
        design_rows.append({"archive_submission_year": year, "platform": platform, "library_strategy": strategy,
                            "library_layout": layout, **counts,
                            "run_weighted_mean_bases_per_archive_count_unit": length_sum / counts["runs_with_count_unit_ratio"] if counts["runs_with_count_unit_ratio"] else None,
                            "aggregate_bases_per_archive_count_unit": counts["complete_record_bases"] / counts["complete_record_reads"] if counts["complete_record_reads"] else None})
    if design_rows:
        with (output / "metadata_designs.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(design_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(design_rows)
    elapsed = time.perf_counter() - start
    report = {
        "schema_version": 1, "input": str(path), "rows_processed": rows, "original_columns": columns,
        "wall_seconds": elapsed, "max_records": max_records,
        "resources": {"parse_and_aggregate_wall_seconds": elapsed,
                      "records_per_second": rows / elapsed if elapsed else None,
                      "input_file_bytes": path.stat().st_size,
                      "peak_rss_mb": process_peak_memory_mb(),
                      "memory_metric": "Peak process working set on Windows; peak RSS on Unix. Includes interpreter/modules and group/organism counters; input records are streamed.",
                      "algorithm": "Single TSV pass; grouped counters, no full table in memory.",
                      "memory_growth": "Proportional to distinct design groups plus organism names; run rows are not retained.",
                      "linear_parse_only_projection_10million_rows_seconds": elapsed / max(1, rows) * 10000000,
                      "projection_caveat": "Extrapolation from observed rows, not a measured10million run; hardware, IO, group cardinality and API acquisition can change scaling."},
        "schema_unavailable_fields": schema_unavailable,
        "missing_counts": dict(missing), "missing_fractions": {k: v / max(rows, 1) for k, v in missing.items()},
        "blank_counts_given_returned_column": {k: v for k, v in missing.items() if k not in schema_unavailable},
        "blank_fractions_given_returned_column": {k: v / max(rows, 1) for k, v in missing.items() if k not in schema_unavailable},
        "missingness_caveat": "missing_counts includes schema absence for backward compatibility. Diagnose archive blank values using blank_counts_given_returned_column; schema_unavailable_fields were not returned/requested and are not archive quality failures.",
        "invalid_counts": dict(invalid), "platform_counts": dict(platforms), "strategy_counts": dict(strategies),
        "layout_counts": dict(layouts), "top_organisms": dict(organisms.most_common(30)),
        "cross_reference_examples": examples, "design_groups": design_rows,
        "scope": "The input file determines archive/query scope. Row count is measured, not assumed to be millions.",
        "limitations": ["Run counts are biased by deposition practices and do not measure design success or independent study counts.",
                        "Year uses public archive submission/release date, not sequencing/collection date.",
                        "base_count/read_count is bases per archive count unit, not necessarily per physical mate. ENA paired read_count can represent spots/pairs (SRR1030394 is an observed example); do not call this mate read length or silently divide all PAIRED rows by2. Raw FASTQ QC supplies actual physical read lengths.",
                        "SRA spots and ENA read_count are stored as distinct fields, but returned read_count semantics may reflect archive spots. No count-unit equivalence is assumed without cross-checking actual records.",
                        "Fields absent from returned/requested TSV columns are schema-unavailable, not evidence of blank/poor archive metadata.",
                        "Streaming groups require memory proportional to groups plus organism names. Full run accession deduplication should be done with SQLite or external sorting before mixing mirrors.",
                        "Library strategy is only a coarse application proxy; variant/expression/assembly labels require study/publication curation."]}
    (output / "metadata_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not quiet:
        print(json.dumps({"rows_processed": rows, "platform_counts": dict(platforms), "output": str(output.resolve()), "resources": report["resources"]}, indent=2))
    return report


def benchmark_metadata(input_path, output_dir, full_report):
    """Measure ordered computational prefixes in fresh processes; reuse full run.

    API/network time is excluded. Prefixes are for performance only, not random
    biological archive samples. Fresh interpreters make peak-RSS comparable.
    """
    output = Path(output_dir)
    records = []
    for count in (1000, 10000, 100000):
        if count >= full_report["rows_processed"]:
            continue
        destination = output / "metadata_scale" / str(count)
        started = time.perf_counter()
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--input", str(Path(input_path).resolve()),
                        "--out", str(destination.resolve()), "--max-records", str(count), "--quiet"], check=True)
        payload = json.loads((destination / "metadata_summary.json").read_text(encoding="utf-8"))
        records.append({"rows": payload["rows_processed"], **payload["resources"],
                        "child_process_wall_seconds_including_startup": time.perf_counter() - started})
    records.append({"rows": full_report["rows_processed"], **full_report["resources"]})
    result = {"schema_version": 1, "input": str(input_path), "benchmarks": records,
              "scope": "Ordered input prefixes in separate processes; computational performance, not probability-sampled biological designs.",
              "caveat": "Parse time excludes network acquisition and result JSON writing. Fresh-process peak RSS includes module startup. Full input measurement is reused from the primary aggregate run; full run need not represent the entire archive."}
    (output / "metadata_scale.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/raw/metadata.tsv")
    parser.add_argument("--out", default="results")
    parser.add_argument("--fetch", action="store_true", help="Acquire ENA metadata first; bounded query, not a probability sample")
    parser.add_argument("--query", default="first_public>=2000-01-01")
    parser.add_argument("--limit", type=int, default=1000000)
    parser.add_argument("--max-records", type=int, help="Ordered computational prefix, intended for scale profiling")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--benchmark", action="store_true", help="Also profile1000/10000/100000-row prefixes in fresh processes")
    arguments = parser.parse_args()
    if arguments.max_records is not None and arguments.max_records < 1:
        parser.error("--max-records must be positive")
    if arguments.benchmark and arguments.max_records is not None:
        parser.error("--benchmark requires the complete acquired input, without --max-records")
    if arguments.fetch:
        fetch_metadata(arguments.input, arguments.query, arguments.limit)
    report = aggregate(arguments.input, arguments.out, arguments.max_records, arguments.quiet)
    if arguments.benchmark:
        benchmark_metadata(arguments.input, arguments.out, report)
