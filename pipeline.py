#!/usr/bin/env python3
"""Streaming QC and reproducible bacterial downsampling experiments.

The dependency-free fallback is an educational seed mapper, NOT minimap2. It
measures conservative exact-anchor-supported reference breadth for ONT. Supply
--sam-illumina/--sam-nanopore from minimap2 to measure standard CIGAR coverage.
No variant, expression, structural-variant or assembly accuracy is inferred.
"""
from __future__ import annotations

import argparse
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys
import time
import tracemalloc

try:
    import numpy as np
except ImportError:
    np = None

ADAPTER = "AGATCGGAAGAGC"
BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
DNA_PATTERN = re.compile(r"[ACGTNacgtn]+")
QUALITY_PATTERN = re.compile(r"[!-~]+")


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="ascii") if path.suffix == ".gz" else path.open(encoding="ascii")


@dataclass(frozen=True)
class Read:
    name: str
    sequence: str
    quality: str


def molecule_id(name):
    """Same key for /1 and /2 and for Illumina whitespace mate fields."""
    return re.sub(r"/[12]$", "", name.split()[0].lstrip("@"))


def fastq_records(path):
    """Strict four-line FASTQ iterator: never accepts silently truncated input."""
    with open_text(path) as handle:
        record = 0
        while True:
            head = handle.readline()
            if not head:
                return
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
            record += 1
            if not seq or not plus or not qual:
                raise ValueError(f"{path}: incomplete FASTQ record {record}")
            head, seq, plus, qual = (s.rstrip("\r\n") for s in (head, seq, plus, qual))
            if not head.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"{path}: malformed FASTQ record {record}")
            if len(seq) != len(qual):
                raise ValueError(f"{path}: sequence/quality length mismatch at record {record}")
            if not seq or DNA_PATTERN.fullmatch(seq) is None:
                raise ValueError(f"{path}: unsupported nucleotide at record {record}")
            if QUALITY_PATTERN.fullmatch(qual) is None:
                raise ValueError(f"{path}: invalid Phred+33 quality at record {record}")
            yield Read(head[1:], seq.upper(), qual)


def read_units(read1, read2=None):
    if read2 is None:
        for read in fastq_records(read1):
            yield (read,)
    else:
        from itertools import zip_longest
        for first, second in zip_longest(fastq_records(read1), fastq_records(read2)):
            if first is None or second is None:
                raise ValueError("Paired FASTQ files have different numbers of records")
            if molecule_id(first.name) != molecule_id(second.name):
                raise ValueError(f"Paired FASTQ IDs disagree: {first.name} / {second.name}")
            yield first, second


def hash_fraction(name, seed):
    digest = hashlib.blake2b(f"{seed}|{molecule_id(name)}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, "big") / 2**64


def selected(name, fraction, seed):
    if not 0 <= fraction <= 1:
        raise ValueError("Downsampling fraction must be in [0,1]")
    return hash_fraction(name, seed) < fraction


def quantile_from_counts(counts, q):
    total = sum(counts.values())
    if not total:
        return None
    target = q * max(0, total - 1)
    seen = 0
    for value, count in sorted(counts.items()):
        seen += count
        if seen > target:
            return value
    return max(counts)


class QC:
    def __init__(self):
        self.reads = self.bases = self.gc = self.n = self.qsum = 0
        self.q20 = self.q30 = self.adapter_reads = 0
        self.expected_errors = 0.0
        self.lengths = Counter()
        self.cycles = [[0, 0] for _ in range(300)]
        self.cycle_sums = np.zeros(300, dtype=np.int64) if np is not None else None
        self.cycle_counts = np.zeros(300, dtype=np.int64) if np is not None else None
        self.sequence_hashes = set()
        self.duplication_sample = 0

    def add(self, read):
        seq, qual = read.sequence, read.quality
        self.reads += 1
        self.bases += len(seq)
        self.gc += seq.count("G") + seq.count("C")
        self.n += seq.count("N")
        self.lengths[len(seq)] += 1
        self.adapter_reads += ADAPTER[:12] in seq
        for char, count in Counter(qual).items():
            q = ord(char) - 33
            self.qsum += q * count
            self.q20 += count if q >= 20 else 0
            self.q30 += count if q >= 30 else 0
            self.expected_errors += 10 ** (-q / 10) * count
        if self.cycle_sums is not None:
            length = min(300, len(qual))
            self.cycle_sums[:length] += np.frombuffer(qual[:300].encode("ascii"), dtype=np.uint8) - 33
            self.cycle_counts[:length] += 1
        else:
            for i, char in enumerate(qual[:300]):
                self.cycles[i][0] += ord(char) - 33
                self.cycles[i][1] += 1
        if self.duplication_sample < 100000:
            self.duplication_sample += 1
            self.sequence_hashes.add(hashlib.blake2b(seq.encode(), digest_size=12).digest())

    def report(self):
        denom = self.bases or 1
        return {
            "reads": self.reads, "bases": self.bases,
            "mean_read_length": self.bases / max(1, self.reads),
            "read_length_min": min(self.lengths, default=0),
            "read_length_p10": quantile_from_counts(self.lengths, .1),
            "read_length_median": quantile_from_counts(self.lengths, .5),
            "read_length_p90": quantile_from_counts(self.lengths, .9),
            "read_length_max": max(self.lengths, default=0),
            "mean_phred": self.qsum / denom,
            "expected_base_error_fraction_from_phred": self.expected_errors / denom,
            "q20_base_fraction": self.q20 / denom, "q30_base_fraction": self.q30 / denom,
            "gc_fraction_acgt": self.gc / max(1, self.bases - self.n),
            "n_base_fraction": self.n / denom,
            "illumina_adapter_prefix12_read_fraction": self.adapter_reads / max(1, self.reads),
            "duplication_fraction_first_100000_sequences": 1 - len(self.sequence_hashes) / max(1, self.duplication_sample),
            "duplication_sample_reads": self.duplication_sample,
            "duplication_caveat": "Identical sequences in a bounded ordered sample; biological and PCR duplicates are not separated.",
            "per_cycle_mean_phred_first_300": [round(float(s) / int(n), 3) if n else None for s, n in
                                                (zip(self.cycle_sums, self.cycle_counts) if self.cycle_sums is not None else self.cycles)],
            "length_histogram": dict(sorted(self.lengths.items())),
        }


def filter_read(read, platform):
    seq, qual = read.sequence, read.quality
    if platform == "illumina":
        adapter_position = seq.find(ADAPTER[:12])
        if adapter_position >= 0:
            seq, qual = seq[:adapter_position], qual[:adapter_position]
        while qual and ord(qual[-1]) - 33 < 20:
            seq, qual = seq[:-1], qual[:-1]
        if len(seq) < 50:
            return None, "length_lt_50_after_adapter_and_q20_tail_trim"
        if seq.count("N") / len(seq) > .02:
            return None, "n_fraction_gt_0.02"
        if sum(qual.encode("ascii")) / len(qual) - 33 < 20:
            return None, "mean_phred_lt_20"
    else:
        # Do not apply short-read tail trimming or adapter sequences to ONT.
        if len(seq) < 500:
            return None, "length_lt_500"
        # Historical ONT Phred calibration varies, so this is a permissive pilot threshold.
        if sum(qual.encode("ascii")) / len(qual) - 33 < 7:
            return None, "mean_phred_lt_7"
    return Read(read.name, seq, qual), None


def load_reference(path):
    names, seqs, chunks = [], [], []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if chunks:
                    seqs.append("".join(chunks).upper())
                names.append(line[1:].split()[0])
                chunks = []
            else:
                chunks.append(line.strip())
    if chunks:
        seqs.append("".join(chunks).upper())
    if not names or len(names) != len(seqs):
        raise ValueError("Reference FASTA must contain named, nonempty contigs")
    offsets, cursor = {}, 0
    for name, seq in zip(names, seqs):
        offsets[name] = cursor
        cursor += len(seq)
    return names, seqs, offsets


def encode_kmer(sequence):
    code = 0
    for base in sequence:
        if base not in BASE_CODE:
            return None
        code = (code << 2) | BASE_CODE[base]
    return code


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def banded_edit_distance(query, target, band):
    """Global Levenshtein distance in O(length*band); returns band+1 if outside."""
    if abs(len(query) - len(target)) > band:
        return band + 1
    inf = band + 1
    previous = {j: j for j in range(min(len(target), band) + 1)}
    for i, a in enumerate(query, 1):
        current = {}
        if i <= band:
            current[0] = i
        for j in range(max(1, i - band), min(len(target), i + band) + 1):
            current[j] = min(previous.get(j, inf) + 1,
                             current.get(j - 1, inf) + 1,
                             previous.get(j - 1, inf) + (a != target[j - 1]))
        previous = current
    return previous.get(len(target), inf)


class SeedMapper:
    """Sampled exact 13-mer index with short-read verification and ONT chains."""
    def __init__(self, seqs, k=13, step=4):
        self.seqs, self.k, self.step = seqs, k, step
        self.offsets = []
        cursor = 0
        packed = []
        for seq in seqs:
            self.offsets.append(cursor)
            for position in range(0, len(seq) - k + 1, step):
                code = encode_kmer(seq[position:position + k])
                if code is not None:
                    packed.append((code << 32) | (cursor + position))
            cursor += len(seq)
        self.length = cursor
        packed.sort()
        self.index = array("Q", packed)
        self.contig_ends = [o + len(s) for o, s in zip(self.offsets, seqs)]

    def hits(self, kmer):
        code = encode_kmer(kmer)
        if code is None:
            return []
        left, right = bisect_left(self.index, code << 32), bisect_left(self.index, (code + 1) << 32)
        # Repetitive seeds do not supply unique placement evidence.
        if right - left > 30:
            return []
        return [self.index[i] & 0xFFFFFFFF for i in range(left, right)]

    def contig_at(self, position):
        return bisect_right(self.offsets, position) - 1

    def map_short(self, sequence):
        candidates = {}
        for strand, query in (("+", sequence), ("-", reverse_complement(sequence))):
            votes = Counter()
            for qpos in range(0, len(query) - self.k + 1, 3):
                for rpos in self.hits(query[qpos:qpos + self.k]):
                    votes[rpos - qpos] += 1
            for start, count in votes.most_common(5):
                if start < 0 or start + len(query) > self.length:
                    continue
                ci = self.contig_at(start)
                if start + len(query) > self.contig_ends[ci]:
                    continue
                local = start - self.offsets[ci]
                target = self.seqs[ci][local:local + len(query)]
                mismatch = sum(a != b for a, b in zip(query, target))
                max_error = max(2, int(len(query) * .15))
                if mismatch <= max_error:
                    distance = mismatch
                else:
                    band = min(max_error, 12)
                    bounded = banded_edit_distance(query, target, band)
                    distance = bounded if bounded <= band else max_error + 1
                if distance <= max_error:
                    key = (start, strand)
                    candidates[key] = (distance, count)
        ordered = sorted((d, -v, pos, strand) for (pos, strand), (d, v) in candidates.items())
        if not ordered:
            return {"mapped": False, "intervals": [], "mapq": None, "error_fraction": None}
        best = ordered[0]
        # Collapse nearby equivalent starts caused by an indel, retaining distant repeats.
        alternative = next((row for row in ordered[1:] if abs(row[2] - best[2]) > 3), None)
        unique = alternative is None or alternative[0] - best[0] >= 3
        return {"mapped": True, "intervals": [(best[2], best[2] + len(sequence))] if unique else [],
                "mapq": None, "unique": unique, "heuristic_score_gap": None if alternative is None else alternative[0] - best[0],
                "strand": best[3], "error_fraction": best[0] / len(sequence)}

    def map_long(self, sequence):
        chains = []
        for strand, query in (("+", sequence), ("-", reverse_complement(sequence))):
            clusters = defaultdict(list)
            for qpos in range(0, len(query) - self.k + 1, 17):
                for rpos in self.hits(query[qpos:qpos + self.k]):
                    # Tolerates accumulated indels through broad diagonal bins.
                    clusters[(self.contig_at(rpos), (rpos - qpos) // 400)].append((qpos, rpos))
            for (ci, diagonal), anchors in clusters.items():
                # Adjacent diagonal bins are combined because an indel can cross a boundary.
                pool = anchors + clusters.get((ci, diagonal + 1), [])
                ordered = sorted(set(pool))
                chain, lastq, lastr = [], -1, -1
                for qpos, rpos in ordered:
                    if qpos > lastq and rpos > lastr:
                        if lastq >= 0 and abs((rpos - lastr) - (qpos - lastq)) > max(60, (qpos - lastq) * .35):
                            continue
                        chain.append((qpos, rpos))
                        lastq, lastr = qpos, rpos
                if len(chain) >= 4 and chain[-1][0] - chain[0][0] >= min(300, len(query) * .3):
                    chains.append((len(chain), chain[-1][0] - chain[0][0], strand, chain))
        if not chains:
            return {"mapped": False, "intervals": [], "mapq": None, "error_fraction": None}
        chains.sort(key=lambda row: (row[0], row[1]), reverse=True)
        best = chains[0]
        rival = next((c for c in chains[1:] if abs(c[3][0][1] - best[3][0][1]) > 1000), None)
        unique = rival is None or best[0] >= rival[0] * 1.5
        # ONLY exact anchor bases contribute. We do not pretend gaps have been aligned.
        intervals = merge_intervals([(rpos, rpos + self.k) for _, rpos in best[3]]) if unique else []
        return {"mapped": True, "intervals": intervals, "mapq": None, "unique": unique,
                "strand": best[2], "anchor_count": best[0],
                "query_anchor_span_fraction": best[1] / len(sequence), "error_fraction": None}


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(merged[-1][1], end)
        else:
            merged.append((start, end))
    return merged


def read_sam(path, offsets, lengths, min_mapq=20):
    """Read primary SAM alignments and retain only M/= /X reference intervals."""
    records = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 11:
                raise ValueError("SAM alignment has fewer than 11 columns")
            name, flag = fields[0], int(fields[1])
            if flag & (256 | 2048):
                continue
            entry = {"mapped": not flag & 4, "intervals": [], "mapq": None,
                     "error_fraction": None, "strand": "-" if flag & 16 else "+"}
            if entry["mapped"]:
                if fields[2] not in offsets:
                    raise ValueError(f"SAM contig {fields[2]} missing from reference")
                mapq = int(fields[4])
                entry["mapq"] = None if mapq == 255 else mapq
                local = int(fields[3]) - 1
                if local < 0:
                    raise ValueError("SAM has non-positive mapped position")
                reference_pos = offsets[fields[2]] + local
                ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
                if not ops or "".join(n + op for n, op in ops) != fields[5]:
                    raise ValueError("Invalid SAM CIGAR")
                aligned_query = 0
                for number, operation in ops:
                    n = int(number)
                    if operation in "M=X":
                        if mapq != 255 and mapq >= min_mapq and not flag & 0x600:
                            entry["intervals"].append((reference_pos, reference_pos + n))
                        reference_pos += n
                        aligned_query += n
                    elif operation in "DN":
                        reference_pos += n
                    elif operation == "I":
                        aligned_query += n
                if reference_pos > offsets[fields[2]] + lengths[fields[2]]:
                    raise ValueError("SAM alignment extends beyond reference contig")
                entry["intervals"] = merge_intervals(entry["intervals"])
                nm = next((int(t[5:]) for t in fields[11:] if t.startswith("NM:i:")), None)
                if nm is not None:
                    entry["error_fraction"] = nm / max(1, aligned_query)
            records[molecule_id(name)].append(entry)
    return records


def coverage_metrics(records, genome_length, gc_windows=None):
    if genome_length <= 0:
        raise ValueError("Reference length must be positive")
    intervals = [interval for record in records for interval in record["intervals"]]
    if any(start < 0 or end < start or end > genome_length for start, end in intervals):
        raise ValueError("Coverage interval outside reference boundaries")
    windows = sorted(gc_windows or [])
    if any(start < 0 or end <= start or end > genome_length or not 0 <= gc <= 1 for start, end, gc in windows):
        raise ValueError("GC window outside reference boundaries or invalid GC fraction")
    if any(windows[i][0] < windows[i-1][1] for i in range(1, len(windows))):
        raise ValueError("GC windows must be disjoint")
    if np is not None:
        delta = np.zeros(genome_length + 1, dtype=np.int32)
        if intervals:
            starts, ends = zip(*intervals)
            np.add.at(delta, np.asarray(starts), 1)
            np.add.at(delta, np.asarray(ends), -1)
        depths = np.cumsum(delta[:-1], dtype=np.int64)
        result = {"coverage_1x": float(np.mean(depths >= 1)), "coverage_5x": float(np.mean(depths >= 5)),
                  "coverage_10x": float(np.mean(depths >= 10)), "mapped_mean_depth_x": float(np.mean(depths))}
        if windows:
            result["gc_dropout"] = [{"start": s, "end": e, "gc_fraction": gc,
                                     "mean_depth": float(np.mean(depths[s:e])),
                                     "uncovered_fraction": float(np.mean(depths[s:e] == 0))} for s, e, gc in windows]
        return result
    # Sparse event sweep uses O(alignment boundaries) memory, not O(genome*reads).
    events = Counter()
    events[0] += 0
    events[genome_length] += 0
    for start, end in intervals:
        events[start] += 1
        events[end] -= 1
    above = Counter()
    covered_bases = current = previous = 0
    gc_depth_sums = [0] * len(windows)
    gc_zero_sums = [0] * len(windows)
    gc_index = 0
    for position, change in sorted(events.items()):
        span = position - previous
        covered_bases += span * current
        for cutoff in (1, 5, 10):
            if current >= cutoff:
                above[cutoff] += span
        while gc_index < len(windows) and windows[gc_index][0] < position:
            window_start, window_end, _ = windows[gc_index]
            overlap = max(0, min(position, window_end) - max(previous, window_start))
            gc_depth_sums[gc_index] += overlap * current
            gc_zero_sums[gc_index] += overlap if current == 0 else 0
            if window_end <= position:
                gc_index += 1
            else:
                break
        current += change
        previous = position
    return {"coverage_1x": above[1] / genome_length, "coverage_5x": above[5] / genome_length,
            "coverage_10x": above[10] / genome_length, "mapped_mean_depth_x": covered_bases / genome_length,
            "gc_dropout": [{"start": start, "end": end, "gc_fraction": gc,
                            "mean_depth": gc_depth_sums[i] / (end - start),
                            "uncovered_fraction": gc_zero_sums[i] / (end - start)}
                           for i, (start, end, gc) in enumerate(windows)]}


def alignment_quality_metrics(records):
    """Read-weighted means of available mapped-primary scores, never truth error."""
    mapqs = [a["mapq"] for a in records if a["mapped"] and a.get("mapq") is not None and a["mapq"] != 255]
    distances = [a["error_fraction"] for a in records if a["mapped"] and a.get("error_fraction") is not None]
    return {"mean_mapq": statistics.mean(mapqs) if mapqs else None,
            "mean_alignment_edit_distance_per_aligned_query_base": statistics.mean(distances) if distances else None,
            "scored_mapq_reads": len(mapqs), "scored_edit_distance_reads": len(distances)}


def platform_metadata(data, platform):
    import csv
    # Broad metadata mining must not accidentally label the pilot with an unrelated
    # accession. A curated manifest of acquired runs takes priority.
    path = data / "selected_runs.tsv"
    if not path.exists():
        path = data / "metadata.tsv"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = row.get("instrument_platform", row.get("platform", "")).upper()
            if (platform == "illumina" and "ILLUMINA" in value) or (platform == "nanopore" and "NANOPORE" in value):
                return {k: row[k] for k in ("run_accession", "study_accession", "sample_accession", "scientific_name", "instrument_platform", "instrument_model", "library_layout") if k in row}
    return {}


def export_comparison(summary, path):
    """Write a scientist-facing two-platform QC/alignment comparison table."""
    import csv
    descriptors = [
        ("Raw physical reads", "reads", lambda d: d["qc"]["reads"], "Both Illumina mates count as separate reads."),
        ("Raw bases examined", "bp", lambda d: d["qc"]["bases"], "Only acquired/examined records, not the entire archived run."),
        ("Raw mean read length", "bp", lambda d: d["qc"]["mean_read_length"], "Measured from physical FASTQ reads; independent of archive count-unit ambiguity."),
        ("Raw median read length", "bp", lambda d: d["qc"]["read_length_median"], "Empirical lower quantile from exact length histogram."),
        ("Raw mean Phred", "Phred", lambda d: d["qc"]["mean_phred"], "Decoded from Phred+33; mean score, not an independently validated sequencing-error rate."),
        ("Raw Q20 base fraction", "fraction", lambda d: d["qc"]["q20_base_fraction"], "Fraction of bases with recorded Phred >=20."),
        ("Raw Q30 base fraction", "fraction", lambda d: d["qc"]["q30_base_fraction"], "Fraction of bases with recorded Phred >=30."),
        ("Raw GC fraction", "fraction", lambda d: d["qc"]["gc_fraction_acgt"], "GC among A/C/G/T bases; N excluded from denominator."),
        ("Raw N fraction", "fraction", lambda d: d["qc"]["n_base_fraction"], "N among all sequenced bases."),
        ("Illumina adapter prefix signal", "fraction of reads", lambda d: d["qc"]["illumina_adapter_prefix12_read_fraction"], "Common Illumina12-base motif; ONT value is a background motif audit, not an ONT adapter screen."),
        ("Identical-sequence sample fraction", "fraction", lambda d: d["qc"]["duplication_fraction_first_100000_sequences"], "First100000 sequences; PCR and biological duplicates cannot be separated."),
        ("Raw identifier units examined", "molecule identifiers", lambda d: d["filtering"]["raw_molecules_examined"], "Both mates share one selection key; a key is not proof of an independent DNA molecule."),
        ("Retained identifier units", "molecule identifiers", lambda d: d["filtering"]["retained_molecules"], "Both mates rejected if either fails eligibility."),
        ("QC-audited retained reads", "reads", lambda d: d["postfilter_qc"]["reads"], "Audit-trimmed sequences describe QC; externally supplied raw SAM is not re-aligned."),
        ("QC-audited mean read length", "bp", lambda d: d["postfilter_qc"]["mean_read_length"], "Internal audit after eligibility/trimming; not external SAM query length."),
        ("QC-audited mean Phred", "Phred", lambda d: d["postfilter_qc"]["mean_phred"], "Decoded from Phred+33; internal retained-sequence audit; ONT is not short-read tail-trimmed."),
        ("Mapped primary reads", "reads", lambda d: d["mapping"]["mapped_reads"], "Includes low-MAPQ mapped primary records."),
        ("Primary mapping rate", "fraction", lambda d: d["mapping"]["mapping_rate"], "Mapped primary reads divided by primary reads of retained identifiers."),
        ("Mean mapped-read MAPQ", "MAPQ", lambda d: d["mapping"]["mean_mapq"], "Available mapped primary MAPQ;255 excluded. Placement confidence, not calibrated truth."),
        ("Mean edit distance per aligned query base", "NM/(M+I+=+X)", lambda d: d["mapping"]["mean_alignment_edit_distance_per_aligned_query_base"], "Sequencing errors and biological reference differences are confounded; not base-error truth."),
        ("Nominal retained depth", "x", lambda d: d["mapping"]["nominal_depth_x"], "Original retained raw bases/reference length in external SAM mode."),
        ("MAPQ-filtered mean mapped depth", "x", lambda d: d["mapping"]["mapped_mean_depth_x"], "Primary matching CIGAR positions after coverage filters; paired overlaps count twice."),
        ("Reference breadth >=1x", "fraction", lambda d: d["mapping"]["coverage_1x"], "Confidence-filtered aligned-read depth endpoint, not variant recall."),
        ("Reference breadth >=5x", "fraction", lambda d: d["mapping"]["coverage_5x"], "At least5 aligned observations; independent-molecule count is unknown."),
        ("Reference breadth >=10x", "fraction", lambda d: d["mapping"]["coverage_10x"], "At least10 aligned observations; no caller or assembly truth benchmark."),
    ]
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "unit", "illumina", "nanopore", "interpretation"])
        for label, unit, getter, interpretation in descriptors:
            writer.writerow([label, unit, getter(summary["datasets"]["illumina"]), getter(summary["datasets"]["nanopore"]), interpretation])


def pilot(args):
    start = time.perf_counter()
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    names, seqs, offsets = load_reference(data / "reference.fasta")
    genome_length = sum(map(len, seqs))
    lengths = dict(zip(names, map(len, seqs)))
    mapper = None
    if not args.sam_illumina or not args.sam_nanopore:
        mapper = SeedMapper(seqs, step=args.index_step)
    reference_gc_windows = []
    for name, seq in zip(names, seqs):
        for pos in range(0, len(seq), 10000):
            window = seq[pos:pos + 10000]
            reference_gc_windows.append((offsets[name] + pos, offsets[name] + pos + len(window),
                                         (window.count("G") + window.count("C")) / len(window)))
    summary = {"schema_version": 1, "reference": {"length_bp": genome_length, "contigs": lengths}, "datasets": {}}
    saturation = {"schema_version": 1, "endpoint": "anchor_supported_reference_breadth",
                  "seeds": [17, 29, 43], "platforms": {}, "datasets": [],
                  "caveat": "Coverage of a bacterial reference is not a variant-calling, expression, assembly or structural-variant accuracy benchmark. Fractions are nested within seed; ranges describe only downsampling randomness, not biological uncertainty."}
    saturation["quality_metric_definitions"] = {
        "mean_mapq": "Unweighted mean available MAPQ over mapped primary reads; MAPQ255 excluded, low MAPQ included. A placement-confidence heuristic, not calibrated truth accuracy.",
        "mean_alignment_edit_distance_per_aligned_query_base": "Unweighted per-read mean NM/(M+I+=+X) for mapped primary SAM records with NM. Sequencing errors and reference biological differences are confounded; deletions enter NM but not query denominator. Educational short fallback uses its verification edit-distance/query-length; long fallback is null.",
        "depth_interpretation": "Per-read scores are reused from a fixed SAM: differences across subset fractions are read-composition/sampling changes, not causal loss of alignment quality from removing unrelated reads."}
    for platform in ("illumina", "nanopore"):
        r1 = data / f"reads_{platform}.fastq.gz"
        r2 = None
        if platform == "illumina" and not r1.exists():
            r1 = data / "reads_illumina_R1.fastq.gz"
            r2 = data / "reads_illumina_R2.fastq.gz"
        if not r1.exists():
            raise FileNotFoundError(f"Missing raw reads: {r1}")
        sam_path = getattr(args, f"sam_{platform}")
        sam = read_sam(sam_path, offsets, lengths, args.min_mapq) if sam_path else None
        limit = args.max_illumina if platform == "illumina" else args.max_nanopore
        qc, filtered_qc, rejection = QC(), QC(), Counter()
        units, raw_units, missing_sam = [], 0, 0
        dataset_start = time.perf_counter()
        for raw in read_units(r1, r2):
            if limit and raw_units >= limit:
                break
            raw_units += 1
            for read in raw:
                qc.add(read)
            filtered = [filter_read(read, platform) for read in raw]
            failures = [reason for _, reason in filtered if reason]
            if failures:
                rejection.update(failures)
                rejection["molecules_rejected"] += 1
                continue
            reads = [read for read, _ in filtered]
            for read in reads:
                filtered_qc.add(read)
            key = molecule_id(reads[0].name)
            if sam is not None:
                alignments = sam.get(key)
                if alignments is None:
                    missing_sam += 1
                    alignments = [{"mapped": False, "intervals": [], "mapq": None, "error_fraction": None} for _ in reads]
            else:
                alignments = [mapper.map_short(r.sequence) if platform == "illumina" else mapper.map_long(r.sequence) for r in reads]
            # SAM must correspond to raw input IDs; preserve aligned-input base count
            # rather than mixing trimmed nominal bases with raw SAM intervals.
            units.append({"id": key, "bases": sum(len(r.sequence) for r in (raw if sam is not None else reads)), "alignments": alignments})
        if missing_sam:
            raise ValueError(f"{platform}: {missing_sam} retained molecule IDs absent from supplied SAM; align the corresponding raw input, not a different subset")
        if not units:
            raise ValueError(f"{platform}: no molecules passed QC")
        if sam is not None:
            method = args.sam_aligner + "_primary_sam_cigar"
            endpoint = f"reference_breadth_primary_alignments_mapq_ge_{args.min_mapq}"
            mapping_caveat = "SAM provenance must specify aligner/version/preset. Primary SAM CIGAR M/= /X bases count toward depth at the MAPQ threshold; MAPQ255 is unavailable and excluded. D/N gaps, secondary/supplementary alignments and preflagged QCFAIL/DUP records do not contribute to coverage; mapped-read counts still include mapped primary QCFAIL/DUP records. No duplicate marking or paired-overlap correction is performed. MAPQ is aligner-specific and not independently calibrated here. NM/(M+I+=+X) is edit distance per aligned query base, includes sequencing errors and biological differences and treats deletions differently from identity. It is not an independent base-error truth benchmark. Internal QC selects accepted molecules from raw input, but trimmed tails/adapters remain in the external raw-SAM intervals. Nominal bases use original supplied input length. Acquired read prefixes can be ordered/biased and are not random samples of the complete archived run."
        else:
            method = "custom_sampled_exact13mer_seeds"
            endpoint = "anchor_supported_reference_breadth" if platform == "nanopore" else "verified_short_read_reference_breadth"
            mapping_caveat = "Educational mapper, not minimap2/BWA. Short reads: exact 13-mer candidates plus mismatch/banded-edit verification; max15% error, top5 candidates/strand, heuristic repeat exclusion. ONT: monotonic exact-anchor chains tolerate indel drift; only exact anchor bases count, so breadth substantially underestimates full alignments. MAPQ is deliberately null because placement probability is not calibrated."
        alignments = [a for u in units for a in u["alignments"]]
        full_coverage = coverage_metrics(alignments, genome_length, reference_gc_windows)
        mapped = sum(a["mapped"] for a in alignments)
        mapqs = [a["mapq"] for a in alignments if a["mapped"] and a["mapq"] is not None]
        errors = [a["error_fraction"] for a in alignments if a["error_fraction"] is not None]
        source = platform_metadata(data, platform)
        summary["datasets"][platform] = {
            "metadata": source, "input_files": [str(r1)] + ([str(r2)] if r2 else []),
            "scope": "real_archived_bacterial_subset", "acquisition_caveat": "Only acquired FASTQ records are analysed; a downloaded prefix may be ordered and cannot represent the complete run without validation.", "max_molecule_limit": limit,
            "qc": qc.report(), "postfilter_qc": filtered_qc.report(),
            "filtering": {"rejection_counts": dict(rejection), "raw_molecules_examined": raw_units,
                          "retained_molecules": len(units), "pair_policy": "Reject both mates if either mate fails" if r2 else "Single file; no paired-end claims",
                          "sam_input_policy": "Accepted molecules from a raw alignment: QC is audited and failing molecules excluded, while SAM intervals and nominal bases refer to original supplied raw input. Trimmed tails/adapters are not re-aligned." if sam is not None else "Custom mapper uses trimmed/filtered sequences"},
            "mapping": {"method": method, "endpoint": endpoint, "caveat": mapping_caveat,
                        "mapped_reads": mapped, "total_reads": len(alignments), "mapping_rate": mapped / len(alignments),
                        "mean_mapq": statistics.mean(mapqs) if mapqs else None,
                        "mean_alignment_edit_distance_per_aligned_query_base": statistics.mean(errors) if errors else None,
                        "nominal_depth_x": sum(u["bases"] for u in units) / genome_length, **full_coverage},
            "processing_seconds": time.perf_counter() - dataset_start,
        }
        points, individual = [], []
        for fraction in (.025, .05, .1, .2, .35, .5, .7, .85, 1.0):
            replicates = []
            for seed in saturation["seeds"]:
                chosen = [u for u in units if selected(u["id"], fraction, seed)]
                rows = [a for u in chosen for a in u["alignments"]]
                metrics = coverage_metrics(rows, genome_length)
                record = {"fraction": fraction, "seed": seed,
                          "nominal_depth_x": sum(u["bases"] for u in chosen) / genome_length,
                          "reads": len(rows), "molecules": len(chosen),
                          "mapping_rate": sum(a["mapped"] for a in rows) / len(rows) if rows else None,
                          **alignment_quality_metrics(rows), **metrics}
                replicates.append(record)
                individual.append(record)
            point = {"fraction": fraction}
            for name in ("nominal_depth_x", "reads", "mapping_rate", "coverage_1x", "coverage_5x", "coverage_10x", "mapped_mean_depth_x",
                         "mean_mapq", "mean_alignment_edit_distance_per_aligned_query_base", "scored_mapq_reads", "scored_edit_distance_reads"):
                values = [r[name] for r in replicates if r[name] is not None]
                point[name + "_mean"] = statistics.mean(values) if values else None
                point[name + "_min"] = min(values) if values else None
                point[name + "_max"] = max(values) if values else None
                point[name + "_replicates_with_value"] = len(values)
            # Expected baseline for ideal independent uniform Poisson sampling.
            depth = point["nominal_depth_x_mean"]
            point["ideal_poisson_breadth_1x"] = 1 - math.exp(-depth)
            points.append(point)
        saturation_depth = None
        for i in range(1, len(points) - 1):
            if (points[i]["coverage_1x_mean"] >= .95 and
                    points[i]["coverage_1x_mean"] - points[i-1]["coverage_1x_mean"] < .01 and
                    points[i+1]["coverage_1x_mean"] - points[i]["coverage_1x_mean"] < .01):
                saturation_depth = points[i]["nominal_depth_x_mean"]
                break
        rule = "At least95% reference breadth>=1x and consecutive absolute gains<1 percentage point; requires two adjacent increments, measured on this subset only."
        saturation["platforms"][platform] = {"metadata": source, "method": method, "endpoint": endpoint,
            "caveat": mapping_caveat, "max_nominal_depth_x": points[-1]["nominal_depth_x_mean"],
            "points": points, "replicates": individual,
            "plateau": {"criterion": rule, "detected": saturation_depth is not None, "depth_x": saturation_depth,
                        "coverage_1x": next((p["coverage_1x_mean"] for p in points if p["nominal_depth_x_mean"] == saturation_depth), None)},
            "mapping_rate_depth_caveat": "Mapping each read uses a fixed reference and fixed algorithm; removing other reads does not causally degrade per-read mapping quality. Changes in sample mapping rate are composition/sampling variability. Breadth and depth are the meaningful depth-dependent endpoints."}
        saturation["datasets"].append({"accession": source.get("run_accession"), "organism": source.get("scientific_name"),
            "platform": platform, "genome_size_mb": genome_length / 1e6, "scope": "bacterial_pilot",
            "points": [{"depth": p["nominal_depth_x_mean"], "mapping_rate": p["mapping_rate_mean"],
                        "breadth_1x": p["coverage_1x_mean"], "breadth_5x": p["coverage_5x_mean"],
                        "breadth_10x": p["coverage_10x_mean"], "mean_mapq": p["mean_mapq_mean"],
                        "mean_alignment_edit_distance_per_aligned_query_base": p["mean_alignment_edit_distance_per_aligned_query_base_mean"]} for p in points],
            "saturation_depth": saturation_depth, "saturation_rule": rule,
            "validation_status": "Real-data reference coverage only; no downstream ground-truth accuracy benchmark"})
    if args.sam_illumina and args.sam_nanopore:
        saturation["endpoint"] = f"reference_breadth_primary_alignments_mapq_ge_{args.min_mapq}"
    # Cross-platform or mixed-method breadth is descriptive, not a superiority test.
    summary["resources"] = {"wall_seconds": time.perf_counter() - start,
                            "coverage_engine": "numpy_difference_arrays" if np is not None else "stdlib_sparse_event_sweep",
                            "python_version": sys.version.split()[0],
                            "peak_rss_mb": process_peak_memory_mb(),
                            "scalability": "FASTQ QC streams; retained per-molecule mappings and sampled reference seed index remain in memory. Full-scale workflow uses minimap2, coordinate-sorted BAM and samtools depth instead."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (out / "saturation.json").write_text(json.dumps(saturation, indent=2), encoding="utf-8")
    export_comparison(summary, out / "qc_alignment_comparison.tsv")
    print(json.dumps({"results": str(out.resolve()), "wall_seconds": summary["resources"]["wall_seconds"],
                      "datasets": {p: {"mapping_rate": v["mapping"]["mapping_rate"], "nominal_depth_x": v["mapping"]["nominal_depth_x"]} for p, v in summary["datasets"].items()}}, indent=2))


def process_peak_memory_mb():
    try:
        import resource
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return value / (1024**2 if sys.platform == "darwin" else 1024)
    except ImportError:
        if sys.platform == "win32":
            import ctypes
            from ctypes import wintypes
            class Counters(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
            counters = Counters()
            counters.cb = ctypes.sizeof(counters)
            get_process = ctypes.windll.kernel32.GetCurrentProcess
            get_process.argtypes = []
            get_process.restype = wintypes.HANDLE
            get_memory = ctypes.windll.psapi.GetProcessMemoryInfo
            get_memory.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            get_memory.restype = wintypes.BOOL
            if get_memory(get_process(), ctypes.byref(counters), counters.cb):
                return counters.PeakWorkingSetSize / 1024**2
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("pilot", help="Real-data QC, alignment audit, nested downsampling and breadth saturation")
    p.add_argument("--data", default="data/raw")
    p.add_argument("--out", default="results")
    p.add_argument("--sam-illumina")
    p.add_argument("--sam-nanopore")
    p.add_argument("--sam-aligner", default="external_SAM", help="Provenance label, e.g. minimap2_2.22_biowasm; never inferred from file extension")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--max-illumina", type=int, default=5000, help="Maximum molecules; 0 processes all acquired reads")
    p.add_argument("--max-nanopore", type=int, default=1000, help="Maximum molecules; 0 processes all acquired reads")
    p.add_argument("--index-step", type=int, default=4)
    args = parser.parse_args()
    if args.command == "pilot":
        pilot(args)


if __name__ == "__main__":
    main()
