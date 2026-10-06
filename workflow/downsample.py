#!/usr/bin/env python3
"""Nested paired molecule hashing, independent of file order and Python hashes."""
import argparse
from contextlib import ExitStack
import gzip
import io
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import molecule_id, read_units, selected


def output_fastq(stack, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    binary = stack.enter_context(gzip.GzipFile(filename=str(path), mode="wb", mtime=0))
    return stack.enter_context(io.TextIOWrapper(binary, encoding="ascii"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--mate")
    parser.add_argument("--output", required=True)
    parser.add_argument("--output-mate")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--fraction", type=float, required=True)
    args = parser.parse_args()
    if bool(args.mate) != bool(args.output_mate):
        parser.error("--mate and --output-mate must be provided together")
    if not 0 <= args.fraction <= 1:
        parser.error("--fraction must be in [0,1]")
    with ExitStack() as stack:
        first = output_fastq(stack, args.output)
        second = output_fastq(stack, args.output_mate) if args.mate else None
        for unit in read_units(args.input, args.mate):
            if selected(molecule_id(unit[0].name), args.fraction, args.seed):
                for read, handle in zip(unit, (first, second)):
                    handle.write(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n")


if __name__ == "__main__":
    main()
