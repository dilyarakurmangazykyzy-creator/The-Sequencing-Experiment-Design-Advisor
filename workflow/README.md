# Standard-tool full-scale workflow

This optional route runs real `fastp`, `minimap2` and `samtools` under Snakemake.
It is separate from the dependency-free pilot and must not be described as run
unless its actual output and version/resource logs exist.

1. Install Mamba/Conda and create an environment from `workflow/environment.yaml`.
2. Put full archived reference/FASTQ files in `data/raw/`, or edit the paths in
   `workflow/config.yaml`. Illumina must supply both mates from the same run.
3. Run from the repository root:

```sh
snakemake --snakefile workflow/Snakefile --cores 4 --use-conda
```

`results/full_workflow_summary.json` contains three independent seeds for six
nested sampling fractions per platform. `alignment_resources.tsv` files report
Snakemake's measured alignment runtime and memory. QC reports specify thresholds.
The configuration fixes the primary-alignment MAPQ filter at 20. `sr` is minimap2's
short-read preset; `map-ont` is the ONT preset. Presets change seed sizes, chaining
and scoring to accommodate platform-specific error modes. Minimizers reduce the
reference/query seed search; seed chains identify candidate placements; dynamic
programming verifies/extends alignments. SAM MAPQ is a heuristic estimate of
placement ambiguity and needs truth-data calibration for a probability claim.

Downsampling hashes normalized molecule IDs with BLAKE2b and an explicit seed;
both mates stay together. An identical read always maps the same way at each
fraction when parameters and reference are fixed. Low depth degrades reference
breadth and downstream evidence, not the intrinsic placement score of a fixed
read. Interpret small mapping-rate changes as sampling/composition variation.

The full route counts aligned M/= /X reference positions, skips deleted/skipped
reference intervals, and excludes secondary/supplementary alignments. Paired
overlaps count twice; PCR duplicates are not identified; no independent variant,
expression or assembly truth set is included. Adding variant recall or assembly
completeness requires a separate endpoint and validation experiment. Conda pins
tool versions but does not lock all transitive package builds; the root pilot
container is the stronger minimal reproducibility option.
